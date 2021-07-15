
**THE SCRIPT TO ANALYZE THE FIRST PILOT's RESULTS***
set more off
clear all

*put your working folder below:
cd C:\Tornado_warnings\Experiment\Alerts_Experiment
log using "./Temp/pilot_analysis.log", replace
set seed 135


*Prepare the treatment characteristics to merge:
import delimited using "./Input/exp_treatments_pilot.csv", varnames(1) clear //I prepare this file separately in R, describes each potential treatment (prior prob+signal information structure)
drop v1

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
rename sequence seq




***CLEANING************************
drop if missing(participant_id) //if no id
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


encode participant_id, gen(subject_id) //as participant_id is initially a string
drop participant_id

save "./Temp/first_pilot_wide.dta", replace

***SANITY CHECKS*******************


**BLIND PROTECTION ANALYSIS (SEPARATELY DUE TO DIFFERENT N of rounds)**
keep subject_id bp_*
reshape long bp_,  i(subject_id) j(round)
rename bp_ bp
gen prob=0.1*round
reg bp prob
xtset subject_id round
xtreg bp prob, fe vce(robust)


**CHANGE TO THE PANEL STRUCTURE FOR THE OTHER TASKS********************
**Panel: participant_id round 
**because all the tasks except the blind protection have the same N of rounds (6) 
use "./Temp/first_pilot_wide.dta", replace

reshape long ip_w_ ip_b_ be_b_ be_w_ wtp_1_ wtp_2_ wtp_3_ wtp_4_ wtp_5_ wtp_6_ wtp_7_ wtp_8_ wtp_9_ wtp_10_ wtp_11_, i(subject_id) j(round)

rename *_ *



***Merging the treatment sequences***
merge m:1 seq round using "./Temp/exp_treatments.dta"
tab _merge
drop _merge

**Calculate the maximum willingness-to-pay for info for each participant/round (wtp):
egen wtp=rowtotal(wtp_1-wtp_11)
replace wtp=0.5*(wtp-1)
sum wtp

**CRUCIAL VARIABLES***

*ip_w - protection in response to seeing a white hint in informed protection (0 - no protection, 1 - protection)
*ip_b - protection in response to seeing a black hint in informed protection (0 - no protection, 1 - protection)
*be_w - belief on prob that the ball is black if seeing a white hint
*be_b - belief on prob that the ball is black if seeing a black hint
*wtp_x - acceptance of price $0.5*(x-1) for information (1 if accepted)
*wtp - revealed wtp for information (in $)


*Saving the cleaned dataset with the panel structure
save "./Output/first_pilot.dta", replace

*Remove the timing information for now
drop *click* *page* history q104-q106 q114 q115




****----RUNNING SOME REGRESSIONS---***********
xtset subject_id round
**informed protection prob response to posterior prob***
xtreg ip_w post_probW, fe vce(robust)
xtreg ip_b post_probB, fe vce(robust)


**belief elicitation reponse to posterior prob (participant FE)**
replace be_w=0.01*be_w
replace be_b=0.01*be_b
xtreg be_w post_probW, fe vce(robust)
xtreg be_b post_probB, fe vce(robust)

*value vs predicted value:
gen value=max(0,-(p*phintWB*loss-(p*phintBB+(1-p)*phintBW)*protectioncost))
sum value
reg wtp value
xtreg wtp value, fe vce(robust)

log close
cmdlog using "./Temp/pilot_analysis.log"

