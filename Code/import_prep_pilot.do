**THE SCRIPT TO IMPORT THE FIRST PILOT'S DATA TO SUBSEQUENTLY MERGE WITH THE MAIN WAVES***
set more off
clear all

*put your working folder below:
cd C:\Tornado_warnings\Experiment\Alerts_Experiment
*log using "./Temp/pilot_analysis.log", replace
set seed 135


*Prepare the treatment characteristics to merge:
import delimited using "./Input/exp_treatments_pilot.csv", varnames(1) clear //I prepare this file separately in R, describes each potential treatment (prior prob+signal information structure)
drop v1
rename snames treatn

save "./Temp/exp_treatments_pilot.dta", replace

import delimited using "./Input/treatment_sequences.csv", varnames(1) clear //I prepare this file separately in Excel, describes treatment order for each group of subjects/session type
rename ïseq seq
reshape long r, i(seq) j(round) //so it is sequence - round structure to make merging with the results file easier
rename r treatn

merge m:1 treatn using "./Temp/exp_treatments_pilot.dta" //merging characteristics of each treatment (prior prob, info structure = N of each gremlin types)
tab _merge
drop _merge
sort seq round

*Create additional variables to describe posterior probabilities:
gen phintWB=(w_gr/tot_gr)
gen phintWW=(tot_gr-bl_gr)/tot_gr
gen phintBB=1-phintWB
gen phintBW=1-phintWW
gen post_probW=p*phintWB/(p*phintWB+(1-p)*phintWW)
gen post_probB=p*phintBB/(p*phintBB+(1-p)*phintBW)

save "./Temp/exp_treatments.dta", replace

*First importing the Qualtrics data
*The results data has variable names incompatible with Stata
*Hence we import the first row of the results file first as data. The first row contains original variable names, which we process to an acceptable format and store in a string 
import delimited using "./Input/Data_FirstPilot.csv", rowrange(1:1) varnames(nonames) clear
local myvarlist
foreach v of varlist _all{
 display `v'
 replace `v'=strtrim(`v')
 replace `v'=subinstr(`v', " (in seconds)", "", .)
 replace `v'=subinstr(`v', " ", "", .)
 replace `v'=substr(`v',3,.)+"_"+substr(`v',1,1) if regexm(`v',"^[0-9]")
 replace `v'=strlower(`v')
 local luck=`v'[1]
 local myvarlist `myvarlist' `luck'
}
*Then we import the data from the results file (ignoring first 3 rows with var names and definition). Then we rename these variables using the previously stored processed names 
import delimited using "./Input/Data_FirstPilot.csv", rowrange(4) varnames(nonames) clear
display "`myvarlist'"
local i=1
foreach v of varlist _all{
 di `i'
 local luck: word `i' of `myvarlist'
 di "`luck'"
 rename `v' `luck'
 local i=`i'+1
}


*Mass renaming some response variables to make reshaping and analysis more convenient:
rename b1_* bp_* //blind protection responses by round
rename be1_b_1* be_b* //beliefs if the gremlin says "black" by round
rename be1_w_15* be_w* //beliefs if the gremlin says "white" by round
rename wtp1_*_1_* wtp_*_* //value responses by round

rename q59_pagesubmit_* bp_time_* //blind protection time to submit
rename q61_pagesubmit_* ip_time_* //informed protection time to submit
rename q79_pagesubmit_* be_time_w_* //belief elicitation time to submit (white hint)
rename q64_pagesubmit_* be_time_b_* //belief elicitation time to submit (black hint)
rename q185_pagesubmit_* wtp_time_* //wtp elicitation time to submit
rename sequence seq




***CLEANING************************
drop if missing(participant_id) //if no id

replace participant_id=participant_id+"PIL"
keep if finished==1 //subjects complete the experiment
drop ipaddress recipient* distributionchannel userlanguage finished //drop unnecessary identifying information
drop round externalreference location* responseid t1r1 sel* treatment_plans ip_wc ip_bc payoffs signal* //other excessive info, relicts of previous experiment edits




*identify correct quiz answers
gen blind_correct=(q157==2)+(q158==3)
tab blind_correct
gen informed_correct=(q120==2)+(q119==2)+(q121==1)+(q118==3)+(q117==4)
tab informed_correct
gen wtpq_correct=(q164==3)+(q166==3)
tab wtpq_correct
gen ncorrect=blind_correct+informed_correct+wtpq_correct
tab ncorrect

*remove participants missing more than 2 questions:
//drop if ncorrect<7
gen pilot=1
save "./Temp/pilot_wide.dta", replace

***SANITY CHECKS*******************


**BLIND PROTECTION PREPATION (SEPARATELY DUE TO DIFFERENT N of rounds)**
* bp - protection decision (0 - do not protect, 1 - protect)

**PREPARE FOR RISK AVERSION CALCULATION CONSISTENT WITH MAIN WAVES
keep participant_id bp_* bp_time_*
reshape long bp_ bp_time_,  i(participant_id) j(round)
rename bp_ bp
rename bp_time_ submittime
gen p=0.1*round
reg bp p

sum submittime

**Identify switchers: change to no protection for higher probabilities
sort participant_id round
gen backswitch=(bp==0)&(bp[_n-1]==1)
replace backswitch=0 if round==1
by participant_id: egen backswitcher=max(backswitch)
by participant_id: egen nbswitches=sum(backswitch)
gen backswitchround=round*backswitch
replace backswitchround=. if backswitchround==0
gen switch=(bp==1)&(bp[_n-1]==0)
replace switch=0 if round==1
by participant_id: egen switcher=max(switch)
gen switchround=round*switch
replace switchround=. if switchround==0

sort participant_id round
list participant_id round bp if backswitcher==1


gen loss=20
gen protectioncost=5
gen bp_val=-((1-bp)*p*loss+bp*protectioncost) //expected costs of blind protection decision
sum bp_val

*if only one back switch and it can be repaired by a single change
*then change the switchround to the 7-total number of round using protection
gen repairable=0
replace repairable=1 if (nbswitches==1)&(bp[_n-1]!=bp)&(bp[_n+1]!=bp)&(round>1)&(round<6)
replace repairable=1 if (nbswitches==1)&(bp[_n-1]!=bp)&(bp==0)&(round==6)



save "./Temp/bp_val_pilot.dta", replace

*Collapsing to have one obs per participant (study participant's characteristics):
collapse (mean) bp submittime (sum) totprot=bp (max) switcher backswitcher nbswitches repairable (min) maxspeed=submittime firstswitch=switchround backswitchround, by(participant_id)

tab firstswitch
replace firstswitch=6-totprot if repairable==1
tab firstswitch

gen backswitcher0=backswitcher
replace backswitcher=0 if repairable==1


tab switcher
tab backswitcher
tab firstswitch //the distribution of the first switching round
tab backswitchround

scatter submittime backswitcher


**Estimate subjects' risk aversion from the blind protection choices:
gen switchprob=0.1*firstswitch
replace switchprob=-0.5 if missing(switchprob)
gen pilot=1
save "./Temp/blind_collapsed_pilot.dta", replace
