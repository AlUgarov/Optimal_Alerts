
**THE SCRIPT TO ANALYZE THE FIRST PILOT's RESULTS***
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
* bp - protection decision (0 - do not protect, 1 - protect)
keep subject_id bp_* bp_time_*
reshape long bp_ bp_time_,  i(subject_id) j(round)
rename bp_ bp
rename bp_time_ submittime
gen p=0.1*round
reg bp p

sum submittime

**Identify switchers: change to no protection for higher probabilities
sort subject_id round
gen switch=(bp==1)&(bp[_n-1]==0)
gen backswitch=(bp==0)&(bp[_n-1]==1)
replace switch=0 if round==1
replace backswitch=0 if round==1
by subject_id: egen backswitcher=max(backswitch)
by subject_id: egen switcher=max(switch)
save "./Temp/blind_prot_temp.dta", replace
gen switchround=round*switch
gen backswitchround=round*backswitch

replace switchround=. if switchround==0
replace backswitchround=. if backswitchround==0

gen loss=20
gen protectioncost=5
gen bp_val=-((1-bp)*p*loss+bp*protectioncost) //expected costs of blind protection decision
sum bp_val

save "./Temp/bp_val.dta", replace

collapse (mean) bp submittime (sum) totprot=bp (max) switcher backswitcher (min) maxspeed=submittime firstswitch=switchround backswitchround, by(subject_id)
tab switcher
tab backswitcher

tab firstswitch //the distribution of the first switching round
tab backswitchround

scatter submittime backswitcher

**Estimate subjects' risk aversion from the blind protection choices:

gen switchprob=0.1*firstswitch

gen byte allprotect=(totprot==5)

save "./Temp/blind_collapsed.dta", replace

use "./Temp/blind_prot_temp.dta", replace

xtset subject_id round
xtreg bp p, fe vce(robust)


**CHANGE TO THE PANEL STRUCTURE FOR THE OTHER TASKS********************
**Panel: participant_id round 
**because all the tasks except the blind protection have the same N of rounds (6) 
use "./Temp/first_pilot_wide.dta", replace


reshape long ip_w_ ip_b_ ip_time_ be_b_ be_w_ be_time_b_ be_time_w_ wtp_1_ wtp_2_ wtp_3_ wtp_4_ wtp_5_ wtp_6_ wtp_7_ wtp_8_ wtp_9_ wtp_10_ wtp_11_ wtp_time_, i(subject_id) j(round)

rename *_ *

rename ip_time time_ip

***Merging the treatment sequences***
merge m:1 seq round using "./Temp/exp_treatments.dta"
tab _merge
drop _merge

**Merging risk aversion vars from the blind protection task
merge m:1 subject_id using "./Temp/blind_collapsed.dta"
drop _merge


**Calculate the maximum willingness-to-pay for info for each participant/round (wtp):
egen wtp=rowtotal(wtp_1-wtp_11)
replace wtp=0.5*(wtp-1)
sum wtp
gen honest_treatment=((bl_gr+w_gr)==0) //only honest gremlins in the group
replace be_w=0.01*be_w //rescale belief elicitation responses to probabilities
replace be_b=0.01*be_b //rescale belief elicitation responses to probabilities


**CRUCIAL VARIABLES***
*ip_ - protection in response to seeing a hint in informed protection (0 - no protection, 1 - protection)
*be_ - belief on prob that the ball is black after seeing a hint in the belief elicitation task
*wtp_x - acceptance of price $0.5*(x-1) for information (1 if accepted)
*wtp - revealed wtp for information (in $)


label var p "Prior prob."
label var phintWB "False neg. rate"
label var phintBW "False pos. rate"

gen ip_val=-(p*(phintWB*(1-ip_w)+phintBB*(1-ip_b))*loss+p*(phintWB*ip_w+phintBB*ip_b)*protectioncost+(1-p)*(phintWW*ip_w+phintBW*ip_b)*protectioncost)
sum ip_val
reg ip_val ncorrect



*Saving the cleaned dataset with the panel structure
save "./Output/first_pilot.dta", replace


*Remove the timing information for now
drop *click* *page* history q104-q106 q114 q115
rename post_probW post_probw
rename post_probB post_probb
reshape long ip_ be_ be_time_ post_prob, i(subject_id round) j(hint) string
gen blackhint=.
replace blackhint=1 if hint=="b"
replace blackhint=0 if hint=="w"
gen question=2*(round-1)+blackhint

sum(be_time_)
sum(time_ip)
label var ip_ "Informed protection"
label var be_ "Belief"
label var blackhint "Gremlin says Black"
label var honest "Honest treatment"
label var post_prob "Posterior prob."


****----RUNNING SOME REGRESSIONS---***********
xtset subject_id question

**--INFORMED PROTECTION--****
**informed protection prob response to posterior prob***
xtreg ip_ post_prob, fe vce(robust)
xtreg ip_ post_prob p blackhint, fe vce(robust)



eststo clear
eststo: xtreg ip_ post_prob, fe vce(robust)
eststo: xtreg ip_ post_prob p blackhint, fe vce(robust)
eststo: xtreg ip_ post_prob if ncorrect>6, fe vce(robust)
eststo: xtreg ip_ post_prob p blackhint if ncorrect>6, fe vce(robust)
esttab using "./Tables/table_ip1.tex", b(%9.3g) t(%9.1f) ar2(%9.2f) label title(Informed Protection) mtitles("All" "All" "Smart" "Smart") star("*" 0.10 "**" 0.05 "***" 0.01) compress nogaps replace

**informed protection prob response to elicited (and posterior) prob***
eststo clear
eststo: xtreg ip_ be_, fe vce(robust)
eststo: xtreg ip_ be_ post_prob, fe vce(robust)
eststo: xtreg ip_ be_ post_prob if ncorrect>6, fe vce(robust)
esttab using "./Tables/table_ip2.tex", b(%9.3g) t(%9.1f) ar2(%9.2f) label title(Informed Protection: Response to Reported Beliefs) mtitles("All" "All" "Smart") star("*" 0.10 "**" 0.05 "***" 0.01) compress nogaps replace

**timing and correlation: see no sizable differences
xtreg ip_ be_, fe vce(robust)
xtreg ip_ be_ if time_ip<20, fe vce(robust)

reg ip_ post_prob, vce(robust)
reg ip_ post_prob if time_ip<20, vce(robust)



**--BELIEF ELICITATION--****
**prepare variables for belief updating responsiveness analysis (Mobius et al, 2011)
gen lt_bel=log(be_/(1-be_))
gen lt_prior=log(p/(1-p))
gen lambda_B=log(phintBB/phintBW)
gen lambda_W=log((1-phintBB)/(1-phintBW))
gen lt_posterior=log(post_prob/(1-post_prob)) //just for verification
gen signalB=blackhint*lambda_B
gen signalW=(1-blackhint)*lambda_W
sum lt_bel lt_posterior lt_prior signalB signalW

reg lt_posterior lt_prior signalB signalW //all the coeffs should be one, zero const

xtreg be_ post_prob, vce(robust)
xtreg be_ post_prob if be_time_>20, vce(robust)

*Testing the accuracy of reported beliefs***
eststo clear
eststo: reg be_ post_prob, vce(robust)
eststo: reg be_ post_prob if ncorrect>6, vce(robust)
eststo: reg be_ post_prob if (honest==0), vce(robust)
esttab using "./Tables/table_be1.tex", b(%9.3g) t(%9.1f) ar2(%9.2f) label title(Belief Elicitation: Belief vs Posterior) mtitles("All" "Not_honest" "Good quiz") star("*" 0.10 "**" 0.05 "***" 0.01) compress nogaps replace


*Testing for biases in reported beliefs***
eststo clear
eststo: reg ip_ post_prob p blackhint, vce(robust)
eststo: xtreg ip_ post_prob p blackhint, fe vce(robust)
eststo: xtreg ip_ post_prob p blackhint if ncorrect>6, fe vce(robust)
esttab using "./Tables/table_be2.tex", b(%9.3g) t(%9.1f) ar2(%9.2f) label title(Belief Elicitation: Determinants) mtitles("OLS" "FE" "Smart, FE") star("*" 0.10 "**" 0.05 "***" 0.01) compress nogaps replace


eststo clear
eststo: reg lt_bel lt_prior signalB signalW, vce(robust)
eststo: xtreg lt_bel lt_prior signalB signalW, fe vce(robust)
eststo: xtreg lt_bel lt_prior signalB signalW if ncorrect>6, fe vce(robust)
esttab using "./Tables/table_be3.tex", b(%9.3g) t(%9.1f) ar2(%9.2f) label title(Belief Elicitation: Decomposition) mtitles("OLS" "FE" "Smart, FE") star("*" 0.10 "**" 0.05 "***" 0.01) compress nogaps replace



**--WTP for Information--***
use "./Output/first_pilot.dta", replace
xtset subject_id round

*theoretical value for a risk-neutral subject:
gen cost_bp=min(p*loss, protectioncost)
gen cost_ip=p*phintWB*loss+(p*phintBB+(1-p)*phintBW)*protectioncost
*gen value=max(0,-(p*phintWB*loss-(p*phintBB+(1-p)*phintBW)*protectioncost))
gen value=max(0,cost_bp-cost_ip)


*belief-adjusted value for a risk-neutral subject:
*gen value_b=max(0,-(p*phintWB*loss-(p*phintBB+(1-p)*phintBW)*protectioncost))



merge m:1 subject_id p using "./Temp/bp_val.dta"

gen info_effect=ip_val-bp_val //expected difference in earnings between informed and blind protection (subject and decision-specific)
sum info_effect

reg info_effect value
reg wtp info_effect
reg wtp info_effect value
reg wtp info_effect value


sum value wtp_time

*Most people answer wtp question in less than 10 seconds
hist wtp_time

reg wtp value
reg wtp value switchprob
reg wtp value totprot
reg wtp value if totprot<5&totprot>0



xtreg wtp value, fe vce(robust)
xtreg wtp value, fe vce(robust)
xtreg wtp value if backswitcher==0, fe vce(robust)
xtreg wtp value if allprotect==0, fe vce(robust)
xtreg wtp value honest_treatment, fe vce(robust)
xtreg wtp value phintWB phintBW, fe vce(robust)
xtreg wtp value honest_treatment phintWB phintBW, fe vce(robust)


eststo clear
eststo: reg wtp value, vce(robust)
eststo: reg wtp value totprot, vce(robust)
eststo: xtreg wtp value, fe vce(robust)
eststo: xtreg wtp value honest_treatment, fe vce(robust)
eststo: xtreg wtp value honest_treatment phintWB phintBW, fe vce(robust)
esttab using "./Tables/table_wtp1.tex", b(%9.3g) t(%9.1f) ar2(%9.2f) label title(WTP for Information) mtitles("OLS" "OLS" "FE" "FE" "FE") star("*" 0.10 "**" 0.05 "***" 0.01) compress nogaps replace

*the speed of answering wtp question does not explain the lack of correlation
reg wtp value
reg wtp value if wtp_time>10


log close
cmdlog using "./Temp/pilot_analysis.log"

