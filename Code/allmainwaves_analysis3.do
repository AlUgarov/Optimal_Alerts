
**ANALYZE THE FIRST WAVE (planned 100 participants)***
set more off
clear all

**Requires: estout, moremata

*!!put your working folder below:
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
gen post_probW=p*phintWB/(p*phintWB+(1-p)*phintWW) //posterior prob that the ball is black if the HINT is white
gen post_probB=p*phintBB/(p*phintBB+(1-p)*phintBW) //posterior prob that the ball is black if the HINT is black

save "./Temp/exp_treatments.dta", replace


*First importing the Qualtrics data
*The results data has variable names incompatible with Stata
*Hence we import the first row of the results file first as data. The first row contains original variable names, which we process to an acceptable format and store in a string 
import delimited using "./Input/Data_All_mainwaves.csv", rowrange(1:1) varnames(nonames) clear
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
import delimited using "./Input/Data_All_mainwaves.csv", rowrange(4:110) varnames(nonames) clear

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

*Checking that the treatment assignment is balanced:
list participant_id seq
tab seq

recode q6 (1=0) (2=1) (4=0), gen(sex)
label var sex "Male"
label define sexes 0 "Non-male" 1 "Male"
label values sex sexes
tab sex

gen age=2021-(2003-75)-q8 //18 means 18 and younger!!

recode q137 (1=1) (2=0) (3=0),gen(stat_educ)
label var stat_educ "Stat. class"
label define stat_educl 0 "No" 1 "Stat. class"
label values stat_educ stat_educl
tab stat_educ

rename q10 educ
label var educ "Education level"
label define educl 1 "No schooling" 2 "Grades 1-12, no diploma" 3 "High school diploma or GED" 4 "Some college, but no degree" 5 "Associate or bachelor's degree" 6 "Graduate or professional degree"
label values educ educl
tab educ


*identify correct quiz answers
gen blind_correct=(q157==2)+(q158==3)+(q189==4)
tab blind_correct
gen informed_correct=(q120==2)+(q119==2)+(q121==1)+(q118==3)+(q117==4)
tab informed_correct
gen add_inf_corect=(q130==1)+(q135==3)+(q136==3)

gen wtpq_correct=(q164==3)+(q166==3)
tab wtpq_correct
gen ncorrect=blind_correct+informed_correct+wtpq_correct
tab ncorrect

gen pilot=0

save "./Temp/mainwaves_wide.dta", replace




**BLIND PROTECTION ANALYSIS**
* bp - protection decision (0 - do not protect, 1 - protect)
keep participant_id bp_* bp_time_* sex age stat_educ ncorrect educ
reshape long bp_ bp_time_,  i(participant_id sex age stat_educ ncorrect) j(round)
rename bp_ bp
rename bp_time_ submittime
gen p=0.05*round
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
save "./Temp/blind_prot_temp.dta", replace


gen loss=20
gen protectioncost=5
gen bp_val=-((1-bp)*p*loss+bp*protectioncost) //expected costs of blind protection decision
sum bp_val

sort participant_id round
list participant_id round bp if backswitcher==1

*if only one back switch and it can be repaired by a single change
*then change the switchround to the 7-total number of round using protection
gen repairable=0
replace repairable=1 if (nbswitches==1)&(bp[_n-1]!=bp)&(bp[_n+1]!=bp)&(round>1)&(round<6)
replace repairable=1 if (nbswitches==1)&(bp[_n-1]!=bp)&(bp==0)&(round==6)


save "./Temp/bp_val.dta", replace




gen college=educ>4



*Collapsing to have one obs per participant (study participant's characteristics):
collapse (mean) bp submittime (first) sex age stat_educ ncorrect college (sum) totprot=bp (max) switcher backswitcher nbswitches repairable (min) maxspeed=submittime firstswitch=switchround backswitchround, by(participant_id)

tab firstswitch
replace firstswitch=. if backswitcher==1
replace firstswitch=7-totprot if repairable==1
tab firstswitch

gen backswitcher0=backswitcher
replace backswitcher=0 if repairable==1

tab switcher
tab backswitcher
tab firstswitch //the distribution of the first switching round
tab backswitchround

gen byte allprotect=(totprot==6)


**Estimate subjects' risk aversion from the blind protection choices:
gen switchprob=0.05*firstswitch
replace switchprob=-0.5 if missing(switchprob)
gen pilot=0





**Merge-in blind protection choices from the pilot
append using "./Temp/blind_collapsed_pilot.dta"

mata:
  function funk(V,p) return(blind_diff(V,30,5,20,p))
  X = st_data(.,("switchprob"))
  Z=J(rows(X),1,0)
  for(i=1; i<=rows(X); i++) {
    p=X[i]
	if (p<0){
	  V=-1
	}
	else {
      mm_root(V=., &funk(), -2, 4, 0.0001, 1000, X[i])
	}
	Z[i]=V
  }
  mata drop funk()
  idx = st_addvar("float", "theta")
  st_store(., "theta", Z)
end

replace theta=. if backswitcher==1
replace theta=9 if (totprot==6)


gen old=age>23&!missing(age) 
label define old_l 0 "$<=$23 yr" 1 "$>$23 yrs"
label values old old_l

gen goodquiz=ncorrect>8
label define goodquiz_l 0 "$>$2 wrong answers" 1 "Good quiz"
label values goodquiz goodquiz_l


//Make the table on distribution of thetas
tab theta
tab pilot



save "./Temp/blind_collapsed.dta", replace

use "./Temp/bp_val.dta", replace
gen pilot=0

*append using "./Temp/bp_val_pilot.dta"
replace pilot=1 if missing(pilot)

keep participant_id p bp bp_val pilot

save "./Temp/bp_val.dta", replace
gen xid=1
collapse (mean) mean=bp (count) n=xid, by(p)
gen se=sqrt(mean*(1-mean)/n)


*Blind protection average response diagram (Fig. 1)
serrbar mean se p, scale (1.96) title("Blind Protection Response") mlwidth(thick) xtitle("Probability of a black ball") ytitle("Proportion of protection choices")
graph export "./Graphs/blind_prot_sta.png", width(1000) height(1000) replace

gen task_type="Blind"
rename p post_prob
keep mean post_prob se task_type
save "./Temp/bp_graph_data.dta", replace




**CHANGE TO THE PANEL STRUCTURE FOR THE OTHER TASKS********************
**Panel: participant_id round 
**because all the tasks except the blind protection have the same N of rounds (6) 
use "./Temp/mainwaves_wide.dta", replace

*add the pilot's data:
*append using "./Temp/pilot_wide.dta"
encode participant_id, gen(subject_id) //as participant_id is initially a string


reshape long ip_w_ ip_b_ ip_time_ be_w_ be_b_ be_time_w_ be_time_b_ wtp_1_ wtp_2_ wtp_3_ wtp_4_ wtp_5_ wtp_6_ wtp_7_ wtp_8_ wtp_9_ wtp_10_ wtp_11_ wtp_time_, i(subject_id participant_id sex age stat_educ) j(round)
rename *_ *
rename ip_time time_ip

***Merging the treatment sequences***
merge m:1 seq round using "./Temp/exp_treatments.dta"
drop _merge


**Merging risk aversion vars from the blind protection task
merge m:1 participant_id using "./Temp/blind_collapsed.dta"
keep if _merge==3
drop _merge

drop if pilot==1 //dropping the pilot

sort participant_id round

**Calculate the maximum willingness-to-pay for info for each participant/round (wtp):
foreach v of varlist wtp_1-wtp_11{
 display `v'
 replace `v'=0 if `v'==-99
}

egen wtp=rowtotal(wtp_1-wtp_11)
replace wtp=0.5*(wtp-1)
replace wtp=0 if wtp<0
sum wtp

gen honest_treatment=((bl_gr+w_gr)==0) //only honest gremlins in the group
replace be_w=0.01*(100-be_w) if pilot==0 //rescale belief elicitation responses to probabilities
replace be_w=0.01*be_w if pilot==1
replace be_b=0.01*be_b //rescale belief elicitation responses to probabilities



**CRUCIAL VARIABLES***
*ip_ - protection in response to seeing a hint in informed protection (0 - no protection, 1 - protection)
*be_ - belief on prob that the ball is black after seeing a hint in the belief elicitation task
*wtp_x - acceptance of price $0.5*(x-1) for information (1 if accepted)
*wtp - revealed wtp for information (in $)
*p - prior probability
*post_prob - posterior probability
*ip_val - expected costs under informed protection given the reported subject's strategy
*ip_val_o - expected costs under informed protection given the cost-minimizing strategy
label var p "Prior prob."
label var phintWB "False neg. rate"
label var phintBW "False pos. rate"

replace ip_w=. if ip_w==-99
replace ip_b=. if ip_b==-99

gen ip_val=-(p*(phintWB*(1-ip_w)+phintBB*(1-ip_b))*loss+p*(phintWB*ip_w+phintBB*ip_b)*protectioncost+(1-p)*(phintWW*ip_w+phintBW*ip_b)*protectioncost)
replace ip_val=. if ip_b==-99
replace ip_val=. if ip_w==-99

sum ip_val



//Calculate the optimal protection strategy based on cost-loss ratio:
gen ip_w_o=0
replace ip_w_o=1 if post_probW>=protectioncost/loss
gen ip_b_o=0
replace ip_b_o=1 if post_probB>=protectioncost/loss


//Calculate exp. costs under the optimal strategy:
gen ip_val_o=-(p*(phintWB*(1-ip_w_o)+phintBB*(1-ip_b_o))*loss+p*(phintWB*ip_w_o+phintBB*ip_b_o)*protectioncost+(1-p)*(phintWW*ip_w_o+phintBW*ip_b_o)*protectioncost)

label var ip_val "Exp. costs"
label var ip_val_o "Optimal exp. costs"

gen ip_val_diff=ip_val-ip_val_o //discrepancy between actual and optimal expected costs of informed protection

replace ip_val_diff=-ip_val_diff

list subject_id p ip_w ip_b ip_val ip_val_o ip_val_diff if abs(ip_val_diff)>3
//large discrepancies emerge when subjects take contrarian actions in the informed protection (e.g. protect when white and do not protect when black)

sum ip_w ip_w_o ip_b ip_b_o if p<0.15&phintWB>0

//Calc exp costs under the optimal strategy for reported(!) beliefs:
gen phintB=p*phintBB+(1-p)*phintBW //ideally we should elicit these probabilities within the experiment
gen phintW=1-phintB

gen ip_w_mu=0
replace ip_w_mu=1 if be_w>=protectioncost/loss
gen ip_b_mu=0
replace ip_b_mu=1 if be_b>=protectioncost/loss

gen ip_val_mu=(loss*(phintW*be_w*(1-ip_w_mu)+phintB*be_b*(1-ip_b_mu))+protectioncost*(phintW*ip_w_mu+phintB*ip_b_mu))
sum ip_val_mu

gen rev_response=(ip_w==1)&(ip_b==0) //subjects protect when white and do not protect when black
tab rev_response


gen treatm_type=0
replace treatm_type=1 if bl_gr>0&w_gr==0
replace treatm_type=2 if w_gr>0&bl_gr==0
replace treatm_type=3 if bl_gr>0&w_gr>0
label define treatm_typel 0 "Honest" 1 "White-eyed only" 2 "Black-eyed only" 3 "All types"
label values treatm_type treatm_typel

*Saving the cleaned dataset with the panel structure
save "./Output/main_waves.dta", replace
use "./Output/main_waves.dta", replace

**Create the summary statistics table:
gen duration_min=duration/60
collapse (first) seq sex age stat_educ college final_payoff duration_min, by(subject_id)

gen seq_type=(seq>3)
tab seq_type seq
gen old=age>23
replace college=1-college

tabstat sex old college stat_educ, by(seq_type) statistics(sum mean) column(statistics) longstub

*collapse (mean) sex old stat_educ college final_payoff duration_min (sum) tsex=sex told=old  tcollege=college tstat_educ=stat_educ
*order tsex sex told old tcollege college tstat_educ stat_educ
*bro


*collapse (mean) sex old stat_educ college final_payoff duration_min (sum) tsex=sex told=old  tcollege=college tstat_educ=stat_educ, by(seq_type)
*order tsex sex told old tcollege college tstat_educ stat_educ
*bro

use "./Output/main_waves.dta", replace

keep subject_id participant_id round time_ip honest ncorrect p phintBW phintWB phintBB rev_response treatm_type
collapse (first) treatm_type participant_id time_ip honest ncorrect p phintBW phintWB phintBB rev_response, by(subject_id round)
label values treatm_type treatm_typel
save "./Temp/base_main_waves.dta", replace

use "./Output/main_waves.dta", replace

*Remove the timing information for now
drop *click* *page* history q104-q106 q114 q115
rename post_probW post_probw
rename post_probB post_probb

keep subject_id participant_id round ip_w ip_b be_w be_b be_time_w be_time_b post_probw post_probb time_ip honest ncorrect p phintBW phintWB phintBB rev_response


*Reshape to (subject_id round hint) long format
reshape long ip_ be_ be_time_ post_prob, i(subject_id round) j(hint) string

*Merging back subject-round variables:
merge m:1 subject_id round using "./Temp/base_main_waves.dta"
keep if _merge==3
drop _merge


**Merging risk aversion vars from the blind protection task
merge m:1 participant_id using "./Temp/blind_collapsed.dta"
keep if _merge==3
drop _merge

gen blackhint=.
replace blackhint=1 if hint=="b"
replace blackhint=0 if hint=="w"
gen question=2*(round-1)+blackhint

gen plevel=round(1000*p) //create integer var to code prior probability

sum(be_time_)
sum(time_ip)
label var ip_ "Informed protection"
label var be_ "Belief"
label var blackhint "Gremlin says Black"
label var honest "Honest treatment"
label var post_prob "Posterior prob."

replace ip_=. if ip_==-99

**Calculate subject-level accuracy of reported beliefs:
gen bel_err=post_prob-be_
label var bel_err "Belief error"

sort subject_id
by subject_id: egen tot_bel_err=sum(abs(bel_err)) //total abs error per subject

sum tot_bel_err, detail
local err_med=r(p50)
hist tot_bel_err
gen accur_bel=tot_bel_err<`err_med' //error is less than the median
label define accur_bel_l 0 "Inaccurate" 1 "Accur. beliefs"
label values accur_bel accur_bel_l

label values stat_educ stat_educl

hist bel_err, title("Errors in elicited beliefs") xtitle("Posterior - Belief") fraction note("By belief elicitation task, no aggregation to round or subjects") color(navy)
graph export "./Graphs/hist_belief_error.png", width(1200) height(800) replace

hist bel_err, title("Errors in elicited beliefs") xtitle("Posterior - Belief") fraction note("Main waves only, excluding certain signals") color(navy)
graph export "./Graphs/hist_belief_error_s3.png", width(1200) height(800) replace

hist bel_err  if pilot==0&abs(0.5-post_prob)<0.499, title("Errors in beliefs, ball color is uncertain") xtitle("Posterior - Belief") fraction note("Main waves only") color(navy)
graph export "./Graphs/hist_belief_error_s4.png", width(1200) height(800) replace


hist bel_err  if pilot==0&abs(0.5-post_prob)>0.499, title("Errors in beliefs, ball color is certain") xtitle("Posterior - Belief") fraction note("Main waves only") color(navy)
graph export "./Graphs/hist_belief_error_s5.png", width(1200) height(800) replace

qui reg be_ post_prob
local r2 : display %5.3f = e(r2)
graph twoway (scatter be_ post_prob, jitter(1)) (lfit be_ post_prob) , title("Belief updating") xtitle("True probability") ytitle("Elicited belief")  note("All obs including pilot") text(0.6 0.5 "R-squared=`r2'", box size(.3cm)) legend(off)
graph export "./Graphs/updating_s1.png", width(1200) height(800) replace

qui reg be_ post_prob if pilot==0&goodquiz==1
local r2 : display %5.3f = e(r2)
graph twoway  (scatter be_ post_prob, jitter(1)) (lfit be_ post_prob) if pilot==0&goodquiz==1, title("Belief updating") xtitle("True probability") ytitle("Elicited belief") legend(off) text(0.6 0.5 "R-squared=`r2'", box size(.3cm)) note("Main waves only, good quiz")
graph export "./Graphs/updating_s2.png", width(1200) height(800) replace

qui reg be_ post_prob if pilot==0&honest==0
local r2 : display %5.3f = e(r2)
graph twoway  (scatter be_ post_prob, jitter(1)) (lfit be_ post_prob) if pilot==0&honest==0, title("Belief updating") xtitle("True probability") ytitle("Elicited belief") legend(off) text(0.6 0.5 "R-squared=`r2'", box size(.3cm)) note("Main waves only, excluding certain signals")
graph export "./Graphs/updating_s3.png", width(1200) height(800) replace

qui reg be_ post_prob if pilot==0&abs(0.5-post_prob)<0.499
local r2 : display %5.3f = e(r2)
graph twoway (scatter be_ post_prob, jitter(1)) (lfit be_ post_prob) if pilot==0&abs(0.5-post_prob)<0.499, title("Belief vs posterior, ball color is uncertain") xtitle("True probability") ytitle("Elicited belief") legend(off) text(0.6 0.5 "R-squared=`r2'", box size(.3cm)) note("Main waves only")
graph export "./Graphs/updating_s4.png", width(1200) height(800) replace

qui reg be_ post_prob if pilot==0&abs(0.5-post_prob)>0.499
local r2 : display %5.3f = e(r2)
graph twoway  (scatter be_ post_prob, jitter(1)) (lfit be_ post_prob) if pilot==0&abs(0.5-post_prob)>0.499, title("Belief vs posterior, ball color is certain") xtitle("True probability") ytitle("Elicited belief") legend(off) text(0.6 0.5 "R-squared=`r2'", box size(.3cm)) note("Main waves only")
graph export "./Graphs/updating_s5.png", width(1200) height(800) replace


**********************************************
****----MAIN REGRESSIONS (beliefs and informed protection)---***********
**********************************************
xtset subject_id question
drop if pilot==1 //now dropping the pilot

save "./Temp/long_ip_dat.dta", replace
use "./Temp/long_ip_dat.dta", replace
gen xid=1
gen post_probi=round(1000*post_prob)

gen blind_prob=0.05*round



*graph twoway (lpoly bp blind_prob, bwidth(0.1) mlwidth(thick)) (lpoly ip post_prob,  mlwidth(thick)), noscatter title("Protection Responses")

*serrbar ip se post_prob if task_type=="Blind", scale (1.96) title("Blind Protection Response") mlwidth(thick) xtitle("Probability of a black ball") ytitle("Proportion of protection choices")

lpoly ip post_prob, bwidth(0.1) ci noscatter title("Informed Protection Response") lineopt(lwidth(1)) mlwidth(thick) xtitle("Posterior probability of a black ball") ytitle("Proportion of protection choices")
graph export "./Graphs/ip_response_lpoly.png", width(1200) height(800) replace

*studying nonsensical responses to certain signals:
by subject_id: egen ip_dev=sd(ip)
gen disagr1=0
gen disagr2=0
replace disagr1=1 if ip_==0&post_prob==1
replace disagr2=1 if ip_==1&post_prob==0

by subject_id: egen totdisagr1=sum(disagr1)

by subject_id: egen totdisagr2=sum(disagr2)

list subject_id totdisagr1 totdisagr2 ip_dev if ip_==0&post_prob==1

gen ip_diff=(ip_-post_prob)^2
bys subject_id: egen tot_diff=sum(ip_diff)
sum tot_diff
list subject_id totdisagr1 totdisagr2 ip_dev tot_diff if ip_==0&post_prob==1

*Notes: nonsensical responses mostly come from subjects choosing different answers to different probabilities, 
* most subjects choosing no protection when prob=1 do not choose protection when prob=0
* seems to be that these subjects are just more random: sum of deviations btw posteriors and responses is much higher for them than for others even excluding that one mistake

collapse (mean) mean=ip_ (count) n=xid, by(post_probi)
gen se=sqrt(mean*(1-mean)/n)

gen task_type="Informed"



gen post_prob=0.001*post_probi
append using "./Temp/bp_graph_data.dta"
serrbar mean se post_prob if task_type=="Informed", scale (1.96) title("Informed Protection Response") mlwidth(thick) xtitle("Posterior probability of a black ball") ytitle("Proportion of protection choices")
graph export "./Graphs/ip_response.png", width(1200) height(800) replace



gen low=mean-1.96*se
gen high=mean+1.96*se
twoway (rcap high low post_prob if task_type=="Blind", lwidth(thick)) (rcap high low post_prob if task_type=="Informed"), title("Protection Response") xtitle("Posterior probability of a black ball") ytitle("Proportion of protection choices") legend(label(1 "Blind") label(2 "Informed"))
graph export "./Graphs/ip_response_comp.png", width(1200) height(800) replace


use "./Temp/long_ip_dat.dta", replace


****-- INFORMED PROTECTION --****

*expected response by prior prob:
gen highprob=p>0.2
label var highprob "p$>$0.2"
label define highprob_l 0 "p $\geq$ 0.2" 1 "p$>$0.2"
label values highprob highprob_l


gen post_prob01=0
replace post_prob01=0 if post_prob<0.2

gen post_prob02=0
replace post_prob02=post_prob if post_prob>=0.2&post_prob<0.4

gen post_prob03=0
replace post_prob03=post_prob if post_prob>=0.4&post_prob<0.6


gen post_prob04=0
replace post_prob02=post_prob if post_prob>=0.6&post_prob<0.8

gen post_prob05=0
replace post_prob05=post_prob if post_prob>=0.8

xtset subject_id
eststo clear
eststo: probit ip_ post_prob0* phintBW phintWB, vce(cluster subject_id)
eststo: probit ip_ post_prob0* phintBW phintWB i.subject_id, vce(cluster subject_id)
eststo: probit ip_ post_prob0* i.highprob##c.phintBW i.highprob##c.phintWB, vce(cluster subject_id)
eststo: probit ip_ post_prob0* i.highprob##c.phintBW i.highprob##c.phintWB i.subject_id, vce(cluster subject_id)
eststo: probit ip_ post_prob0* i.blackhint##c.phintBW i.blackhint##c.phintWB, vce(cluster subject_id)
eststo: probit ip_ post_prob0* i.blackhint##c.phintBW i.blackhint##c.phintWB i.subject_id, vce(cluster subject_id)
esttab using "./Tables/table_ip5.tex", b(%9.3g) se(%9.1f) ar2(%9.2f) label indicate(Subject FE=*.subject_id) drop(post_prob0*) mtitles("" "" "" "" "" "") title(Informed protection by prior) star("*" 0.10 "**" 0.05 "***" 0.01) nobaselevels compress nogaps replace




gen post_prob2=post_prob^2

eststo clear
eststo: probit ip_ post_prob post_prob2 phintBW phintWB , vce(cluster subject_id)
eststo: probit ip_ post_prob post_prob2 phintBW phintWB if plevel==100, vce(cluster subject_id)
eststo: probit ip_ post_prob post_prob2 phintBW phintWB if plevel==200, vce(cluster subject_id)
eststo: probit ip_ post_prob post_prob2 phintBW phintWB if plevel==300, vce(cluster subject_id)
eststo: probit ip_ post_prob post_prob2 phintBW phintWB if plevel==500, vce(cluster subject_id)
esttab using "./Tables/table_ip4.tex", b(%9.3g) se(%9.1f) ar2(%9.2f) label title(Informed protection by prior) mtitles("All" "0.1" "0.2" "0.3" "0.5") star("*" 0.10 "**" 0.05 "***" 0.01) nobaselevels compress nogaps replace


gen highprobBW=highprob*phintBW
gen highprobWB=highprob*phintWB

gen blackhintBW=blackhint*phintBW
gen blackhintWB=blackhint*phintWB
gen statBW=stat_educ*phintBW
gen statWB=stat_educ*phintWB


eststo clear
eststo: semipar ip_ phintBW phintWB, nonpar(post_prob)
eststo: semipar ip_ highprob phintBW highprobBW phintWB highprobWB, nonpar(post_prob) 
eststo: semipar ip_ blackhint phintBW phintWB blackhintBW blackhintWB, nonpar(post_prob)
eststo: semipar ip_ stat_educ phintBW phintWB statBW statWB, nonpar(post_prob)
esttab using "./Tables/table_ip5_semi.tex", b(%9.3g) se(%9.1f) ar2(%9.2f) label mtitles("" "" "" "") title(Informed protection by prior) star("*" 0.10 "**" 0.05 "***" 0.01) nobaselevels compress nogaps replace




****-- BELIEF ELICITATION --****
*Testing the accuracy of reported beliefs***
eststo clear
eststo: reg bel_err phintWB phintBW, vce(cluster subject_id)
eststo: reg bel_err i.plevel phintWB phintBW, vce(cluster subject_id)
eststo: reg bel_err i.goodquiz##c.phintWB i.goodquiz##c.phintBW, vce(cluster subject_id)
eststo: reg bel_err i.goodquiz##i.plevel i.goodquiz##c.phintWB i.goodquiz##c.phintBW, vce(cluster subject_id)
eststo: reg bel_err i.stat_educ##c.phintWB i.stat_educ##c.phintBW, vce(cluster subject_id)
eststo: reg bel_err i.stat_educ##i.plevel i.stat_educ##c.phintWB i.stat_educ##c.phintBW, vce(cluster subject_id)
esttab using "./Tables/table_be2.tex", b(%9.3g) se(%9.1f) ar2(%9.2f) label title(Belief Elicitation: Discrepancy) mtitles("" "" "" "" "" "") star("*" 0.10 "**" 0.05 "***" 0.01) indicate(Prior prob dummies = *.plevel) nobaselevels compress nogaps replace


**prepare variables for belief updating responsiveness analysis (Mobius et al, 2011)
gen lt_bel=log(be_/(1-be_))
gen lt_prior=log(p/(1-p))
gen lambda_B=log(phintBB/phintBW)
gen lambda_W=log((1-phintBB)/(1-phintBW))
gen lt_posterior=log(post_prob/(1-post_prob)) //just for verification
gen signalB=blackhint*lambda_B
gen signalW=(1-blackhint)*lambda_W
gen signal=signalB+signalW
sum lt_bel lt_posterior lt_prior signalB signalW

label var lt_prior "Prior"
label var signal "Signal"
*Decomposition of belief updating (coeffs should be all ones)***
eststo clear
eststo: reg lt_bel lt_prior signal, noconstant vce(cluster subject_id)
eststo: xtreg lt_bel lt_prior signal, fe vce(cluster subject_id)
eststo: reg lt_bel lt_prior signal i.goodquiz#c.lt_prior i.goodquiz#c.signal, noconstant vce(cluster subject_id)
eststo: xtreg lt_bel lt_prior signal i.goodquiz#c.lt_prior i.goodquiz#c.signal, fe vce(cluster subject_id)
eststo: reg lt_bel lt_prior signal i.stat_educ#c.lt_prior i.stat_educ#c.signal, noconstant vce(cluster subject_id)
eststo: xtreg lt_bel lt_prior signal i.stat_educ#c.lt_prior i.stat_educ#c.signal, fe vce(cluster subject_id)
esttab using "./Tables/table_be3.tex", b(%9.3g) t(%9.1f) ar2(%9.2f) label drop(_cons) mtitles("OLS" "FE" "OLS" "FE" "OLS" "FE") star("*" 0.10 "**" 0.05 "***" 0.01) note("Decomposition works only for imperfect signals") nobaselevels compress nogaps replace

gen psignalblack=p*phintBB+(1-p)*phintBW
gen psignalwhite=1-psignalblack
gen pifsignblack=p*phintBB/psignalblack


**********************************************
****-- Main Regressions: WTP FOR INFORMATION --****
**********************************************
use "./Output/main_waves.dta", replace
xtset subject_id round

merge m:1 subject_id using "./Temp/bel_accuracy.dta" //average belief accuracy by subject
drop _merge


merge m:1 participant_id p using "./Temp/bp_val.dta" //risk aversion, blind prot choices and demographic vars
drop if _merge==2
drop _merge

drop if pilot==1 //dropping the pilot

*calculate theoretical value of signal for a risk-neutral subject:
gen cost_bp=min(p*loss, protectioncost) //expected blind protection cost
gen false_pos=(1-p)*phintBW*protectioncost //expected protection costs for false positives
gen true_pos=p*phintBB*protectioncost //expected protection costs for true positives

gen false_neg=p*phintWB*loss //expected false positive costs
gen cost_ip=false_neg+false_pos+true_pos //total informed protection costs
gen pos_cost=false_pos+true_pos //total protection costs in response to positive signals
gen value=max(0, cost_bp-cost_ip) //theoretical value for a risk-neutral subject

gen freqBW=(1-p)*phintBW
gen freqBB=p*phintBB
gen freqWB=p*phintWB


label var false_pos "FP costs"
label var false_neg "FN costs"
label var true_pos "TP costs"
label var pos_cost "Pos. signal costs"
label var cost_bp "BP costs"
label var phintWB "FN rate"
label var phintBW "FP rate"




*uniform CRRA values to calculate the signal's values later
gen theta1=0.5
gen theta2=1
gen theta3=1.5
gen theta4=2.5

*calculate theoretical value for a risk-averse subject (both heterogeneous and uniform CRRA coeffs):
qui mata:
  //mata drop val_cara()
  //mata drop myfunc2()
  function myfunc2(V,p,pWW,pBB,thet) return(infoval_diff(V,30,5,20,p,pWW,pBB,thet))
  function val_cara(V,p,pWW,pBB,thet) return(infoval_diff_cara(V,30,5,20,p,pWW,pBB,thet))
  X = st_data(.,("p","phintWW","phintBB","theta"))
  Z=J(rows(X),1,0)
  
  Z1=J(rows(X),1,0)
  Z2=J(rows(X),1,0)
  Z3=J(rows(X),1,0)
  Z4=J(rows(X),1,0)
  Z5=J(rows(X),1,0)
  Z6=J(rows(X),1,0)
  
  for(i=1; i<=rows(X); i++) {
    p=X[i,1]
	phintWW=X[i,2]
	phintBB=X[i,3]
	theta=X[i,4]
	if ((theta<-0.99)||(theta==4)){
	  V=-99
	}
	else {
      mm_root(V=., &myfunc2(), -2, 5, 0.0001, 1000, p,phintWW,phintBB,theta)
	}
	Z[i]=V
	mm_root(V1=., &myfunc2(), -2, 5, 0.0001, 1000, p,phintWW,phintBB,0.5)
	mm_root(V2=., &myfunc2(), -2, 5, 0.0001, 1000, p,phintWW,phintBB,1)
	mm_root(V3=., &myfunc2(), -2, 5, 0.0001, 1000, p,phintWW,phintBB,1.5)
	mm_root(V4=., &myfunc2(), -2, 5, 0.0001, 1000, p,phintWW,phintBB,2.5)
	
	mm_root(V5=., &val_cara(), -2, 5, 0.0001, 1000, p,phintWW,phintBB,0.05)
	mm_root(V6=., &val_cara(), -2, 5, 0.0001, 1000, p,phintWW,phintBB,0.3)
	Z1[i]=V1
	Z2[i]=V2
	Z3[i]=V3
	Z4[i]=V4
	Z5[i]=V5
	Z6[i]=V6
  }
  mata drop myfunc2()
  idx = st_addvar("float", "value_ra")
  st_store(., "value_ra", Z)
  
  idx = st_addvar("float", "value_ra1")
  st_store(., "value_ra1", Z1)
  
  idx = st_addvar("float", "value_ra2")
  st_store(., "value_ra2", Z2)
  
  idx = st_addvar("float", "value_ra3")
  st_store(., "value_ra3", Z3)
  
  idx = st_addvar("float", "value_ra4")
  st_store(., "value_ra4", Z4)
  
  idx = st_addvar("float", "value_cara1")
  st_store(., "value_cara1", Z5)
  
  idx = st_addvar("float", "value_cara2")
  st_store(., "value_cara2", Z6)
  
end

replace value_ra=. if ((backswitcher==1))
replace value_ra=0 if value_ra<0
replace value_ra1=0 if value_ra1<0
replace value_ra2=0 if value_ra2<0
replace value_ra3=0 if value_ra3<0
replace value_ra4=0 if value_ra4<0
replace value_cara1=0 if value_cara1<0
replace value_cara2=0 if value_cara2<0
replace value_ra=. if theta==-1


gen wtp_diffra=wtp-value_ra
replace wtp_diffra=0 if wtp_diffra<0



*risk-averse indicator:
gen risk_averse=theta>0&backswitcher==0
replace risk_averse=. if backswitcher==1

label var risk_averse "Risk-averse indicator"
label define risk_averse_l 0 "Not risk averse" 1 "Risk-averse"
label values risk_averse risk_averse_l


gen phintWBs=phintWB*loss
gen phintBWs=phintBW*protectioncost
label var phintWBs "False-neg. prob. x Loss"
label var phintBWs "False-neg. prob. x Prot. cost"


replace bp_val=-bp_val

gen info_effect=bp_val+ip_val //expected difference in earnings between informed and blind protection (subject and decision-specific)
replace info_effect=0 if info_effect<0 //because wtp is never negative


gen value_mu=max(0, bp_val-ip_val_mu) //expected diff in earnings between IP and BP accounting for actual subjects' beliefs

*Calculate discrepancies between WTP and theoretical value for different risk aversion levels:
gen wtp_diff=wtp-value

label values accur_bel accur_bel_l //restoring lost value labels



*Histograms:
hist ip_val_diff, title("Distribution of expected costs discrepancies") xtitle("Discrepancy") fraction note("Difference between optimal and actual expected costs in the IP treatment (by round)") color(navy)
graph export "./Graphs/hist_costs_discr.png", width(1200) height(800) replace

hist wtp_diff, title("Distribution of WTP discrepancies (WTP - Value)") xtitle("USD") fraction note("Difference between stated wtp and theoretical value for a risk-neutral subject (each choice=obs)") color(navy)
graph export "./Graphs/hist_WTP_discr1.png", width(1200) height(800) replace

sort subject_id
gen wtp_diff_abs=abs(wtp_diff)
by subject_id: egen totwtp_diff=sum(wtp_diff_abs)
replace totwtp_diff=(1/6)*totwtp_diff

hist totwtp_diff, title("Distribution of WTP discrepancies (WTP - Value)") xtitle("USD") fraction note("Average absolute deviation between stated wtp and theoretical value for a risk-neutral subject, by subject") color(navy)
graph export "./Graphs/hist_WTP_discr2.png", width(1200) height(800) replace


scatter wtp value, title("Theoretical vs actual WTP") xtitle("Theoretical WTP") ytitle("Actual WTP") jitter(2) 
graph export "./Graphs/WTP_value_scatter.png", width(1200) height(800) replace

heatplot wtp value, normalize backfill colors(HCL blues, reverse) title("Theoretical vs actual WTP") xtitle("Theoretical WTP") ytitle("Actual WTP") levels(10) bins(10) keylabels(minmax)
graph export "./Graphs/WTP_value_heat.png", width(1200) height(800) replace

heatplot wtp value_ra, normalize backfill colors(HCL blues, reverse) title("Theoretical with risk aversion vs actual WTP") xtitle("Theoretical WTP") ytitle("Actual WTP") levels(10) bins(10) keylabels(minmax)
graph export "./Graphs/WTP_value_heat_ra.png", width(1200) height(800) replace

*Note: correlation wtp value - 0.27; correlation wtp value_ra - 0.24

heatplot wtp value, colors(s2) title("Theoretical vs actual WTP") xtitle("Theoretical WTP") ytitle("Actual WTP") discrete levels(1) scatter sizeprop keylabels(none)
graph export "./Graphs/WTP_value_heat_scatter.png", width(1200) height(800) replace

gen highprob=p>0.2
label var highprob "p$>$0.2"
label define highprob_l 0 "p $\geq$ 0.2" 1 "p$>$0.2"
label values highprob highprob_l

gen plevel=round(1000*p)
gen risk_loving=theta<-0.001
replace risk_loving=. if missing(theta)
label var risk_loving "Risk loving"
label define risklov_l 0 "Not risk-loving" 1 "Risk loving"
label value risk_loving risklov_l

gen risk_missing=missing(theta)
label var risk_missing "Risk loving"
label define riskmiss_l 0 "RA measured" 1 "No risk av. measure"
label value risk_missing riskmiss_l

replace risk_loving=0 if risk_missing==1
replace risk_averse=0 if risk_missing==1

gen risk_neutral=(risk_averse+risk_loving==0)&!missing(theta)
gen model_subject=risk_neutral+accur_bel==2

label define model_subject_l 0 "Not" 1 "Risk-neutral and accur."
label value model_subject model_subject_l

gen nmodel_subject=1-model_subject
*replace nmodel_subject=. if missing(risk_loving)
label define nmodel_subject_l 0 "Risk-neutral and accur." 1 "Not risk-neutral and accurate"
label value nmodel_subject nmodel_subject_l


gen risk_pref=0 if risk_neutral==1

replace risk_pref=1 if risk_loving==1
replace risk_pref=2 if risk_averse==1
replace risk_pref=3 if missing(theta)
label var risk_pref "Risk preferences"
label define risk_prefl 0 "Risk-neutral" 1 "Risk-loving" 2 "Risk-averse" 3 "No risk av. measure"
label value risk_pref risk_prefl

gen accur_bel1=1-accur_bel
label var accur_bel1 "Beliefs accuracy"
label def accur_bel1l 0 "Accur. beliefs" 1 "Inaccurate beliefs"
label value accur_bel1 accur_bel1l


**Baseline WTP difference: risk-aversion and belief accuracy:
eststo clear
eststo: reg wtp_diff false_pos false_neg, vce(cluster subject_id)
eststo: reg wtp_diff i.risk_averse##c.false_pos i.risk_averse##c.false_neg i.risk_loving##c.false_pos i.risk_loving##c.false_neg i.risk_missing##c.false_pos i.risk_missing##c.false_neg , vce(cluster subject_id)
eststo: reg wtp_diff i.accur_bel##c.false_pos i.accur_bel##c.false_neg, vce(cluster subject_id)
eststo: reg wtp_diff i.accur_bel1#i.risk_pref##c.false_pos i.accur_bel1#i.risk_pref##c.false_neg, vce(cluster subject_id)
esttab using "./Tables/table_wtpdiff_01rfull.tex", b(%9.3g) se(%9.1f) ar2(%9.2f) label title(WTP for Information (Discrepancy)) mtitles("" "" "" "") star("*" 0.10 "**" 0.05 "***" 0.01) nobaselevels compress nogaps replace

esttab using "./Tables/table_wtpdiff_01r.tex", b(%9.3g) se(%9.1f) ar2(%9.2f) label title(WTP for Information (Discrepancy)) mtitles("" "" "" "") indicate("Risk aversion#Belief accuraccy dummies = *.accur_bel1#*.risk_pref") star("*" 0.10 "**" 0.05 "***" 0.01) nobaselevels compress nogaps replace



*Demographic variables:
eststo clear
eststo: reg wtp_diff false_pos false_neg, vce(cluster subject_id)
eststo: reg wtp_diff i.sex##c.false_pos i.sex##c.false_neg, vce(cluster subject_id)
eststo: reg wtp_diff i.sex##i.plevel i.sex##c.false_pos i.sex##c.false_neg, vce(cluster subject_id)
eststo: reg wtp_diff i.stat_educ##c.false_pos i.stat_educ##c.false_neg, vce(cluster subject_id)
eststo: reg wtp_diff i.stat_educ##i.plevel i.stat_educ##c.false_pos i.stat_educ##c.false_neg, vce(cluster subject_id)
eststo: reg wtp_diff i.old##c.false_pos i.old##c.false_neg, vce(cluster subject_id)
eststo: reg wtp_diff i.old##i.plevel i.old##c.false_pos i.old##c.false_neg, vce(cluster subject_id)
eststo: reg wtp_diff i.goodquiz##c.false_pos i.goodquiz##c.false_neg, vce(cluster subject_id)
eststo: reg wtp_diff i.goodquiz##i.plevel i.goodquiz##c.false_pos i.goodquiz##c.false_neg, vce(cluster subject_id)
esttab using "./Tables/table_wtpdiff_02.tex", b(%9.3g) se(%9.1f) ar2(%9.2f) label indicate(Prior dummies=*.plevel) title(WTP for Information (Discrepancy, demographic variables)) mtitles("" "" "" "" "" "" "" "" "") star("*" 0.10 "**" 0.05 "***" 0.01) nobaselevels compress nogaps replace


*By prior prob:
eststo clear
eststo: reg wtp_diff false_pos false_neg if plevel==100, vce(cluster subject_id)
eststo: reg wtp_diff false_pos false_neg if plevel==200, vce(cluster subject_id)
eststo: reg wtp_diff false_pos false_neg if plevel==300, vce(cluster subject_id)
eststo: reg wtp_diff false_pos false_neg if plevel==500, vce(cluster subject_id)
esttab using "./Tables/table_wtpdiff_03.tex", b(%9.3g) se(%9.1f) ar2(%9.2f) label title(WTP for Information (Discrepancy, by prior)) mtitles("0.1" "0.2" "0.3" "0.5") star("*" 0.10 "**" 0.05 "***" 0.01) nobaselevels compress nogaps replace


*Testing heterogeneity
reg wtp_diff i.plevel##c.false_pos i.plevel##c.false_neg, vce(cluster subject_id)
contrast plevel plevel#c.false_pos plevel#c.false_neg, overall


*By prior prob 2:
eststo clear
eststo: reg wtp_diff phintBW phintWB, vce(cluster subject_id)
eststo: reg wtp_diff phintBW phintWB if plevel==100, vce(cluster subject_id)
eststo: reg wtp_diff phintBW phintWB if plevel==200, vce(cluster subject_id)
eststo: reg wtp_diff phintBW phintWB if plevel==300, vce(cluster subject_id)
eststo: reg wtp_diff phintBW phintWB if plevel==500, vce(cluster subject_id)
esttab using "./Tables/table_wtpdiff_06.tex", b(%9.3g) se(%9.1f) ar2(%9.2f) label title(WTP for Information (Discrepancy, by prior)) mtitles("All" "0.1" "0.2" "0.3" "0.5") star("*" 0.10 "**" 0.05 "***" 0.01) nobaselevels compress nogaps replace



*By prior prob:
eststo clear
eststo: tobit wtp highprob phintBW i.highprob#c.phintBW phintWB i.highprob#c.phintWB, ll(0) ul(5)
eststo: tobit wtp highprob phintBW i.highprob#c.phintBW phintWB i.highprob#c.phintWB if goodquiz==1, ll(0) ul(5)
eststo: tobit wtp highprob phintBW i.highprob#c.phintBW phintWB i.highprob#c.phintWB if stat_educ==1, ll(0) ul(5)
*eststo: tobit value highprob phintBW i.highprob#c.phintBW phintWB i.highprob#c.phintWB, ll(0) ul(5)
esttab using "./Tables/table_wtp_val0.tex", b(%9.3g) se(%9.1f) ar2(%9.2f) label title() mtitles("All" "Good quiz only" "Stat. classes") star("*" 0.10 "**" 0.05 "***" 0.01) nobaselevels compress nogaps replace


*By prior prob:
eststo clear
eststo: tobit wtp highprob phintBW i.highprob#c.phintBW phintWB i.highprob#c.phintWB, ll(0) ul(5)
eststo: tobit wtp i.goodquiz##c.highprob i.goodquiz##c.phintBW i.goodquiz##i.highprob#c.phintBW i.goodquiz##c.phintWB i.goodquiz##i.highprob#c.phintWB, ll(0) ul(5)
eststo: tobit wtp i.stat_educ##c.highprob i.stat_educ##c.phintBW i.stat_educ##i.highprob#c.phintBW i.stat_educ##c.phintWB i.stat_educ##i.highprob#c.phintWB, ll(0) ul(5)
*eststo: tobit value highprob phintBW i.highprob#c.phintBW phintWB i.highprob#c.phintWB, ll(0) ul(5)
esttab using "./Tables/table_wtp_val0ext.tex", b(%9.3g) se(%9.1f) ar2(%9.2f) label title() mtitles("" "" "") star("*" 0.10 "**" 0.05 "***" 0.01) nobaselevels compress nogaps replace


*Create variables for the next table:
gen ps2=p^2
label var ps2 "p squared"
label var freqBW "FP total prob."
label var freqWB "FN total prob."
gen freqBWs=freqBW^2
gen freqWBs=freqWB^2
label var freqBWs "FP total prob. sq."
label var freqWBs "FN total prob. sq."

*Interaction stuff:
eststo clear
eststo: tobit wtp p freqBW freqWB phintBW phintWB, ll(0) ul(5)
eststo: tobit wtp p freqBW freqWB phintBW phintWB ps2 freqBWs freqWBs, ll(0) ul(5)
eststo: tobit wtp p freqBW freqWB i.stat_educ##c.phintBW i.stat_educ##c.phintWB, ll(0) ul(5)
eststo: tobit wtp p freqBW freqWB i.stat_educ##c.phintBW i.stat_educ##c.phintWB ps2 freqBWs freqWBs, ll(0) ul(5)
*eststo: tobit value p freqBW freqWB phintBW phintWB, ll(0) ul(5)
*eststo: tobit value p freqBW freqWB phintBW phintWB ps2 freqBWs freqWBs, ll(0) ul(5)
esttab using "./Tables/table_wtp_val1.tex", b(%9.3g) se(%9.1f) ar2(%9.2f) label drop(p freqBW freqWB ps2 freqWBs) indicate(With squares=freqBWs) title() note("Controlling for priors and total probabilities of false-posiive and false-negative outcomes. Standard errors in parentheses.") mtitles("" "" "" "") star("*" 0.10 "**" 0.05 "***" 0.01) nobaselevels compress nogaps replace


*Paying positive amounts for signals not affecting their protection choices:
gen pay_notuse=(wtp>0)&(ip_b==ip_w)
gen use_signal=(ip_b>ip_w)
tab pay_notuse //many choices are inconsistent

//most subjects make at least one inconsistent choice here:
sort subject_id
by subject_id: egen totpay_notuse=sum(pay_notuse)
tab totpay_notuse

*sign-switching patterns: high sensitivity to FP when the prior prob is low and high sensit to FN when the prior is high
bys plevel: reg wtp_diff false_pos false_neg, vce(cluster subject_id)

*checking if the sign-switching is due to treatments in which it is better not to respond to signals
bys plevel: reg wtp_diff false_pos false_neg if value>0, vce(cluster subject_id)


*measuring sensitivity of false-negative rates:
*gen highprob=p>0.2
collapse (mean) wtp_diff (sd) wtp_diffsd=wtp_diff (count) nobs=wtp_diff, by(highprob phintBW phintWB)
gen sign=(wtp_diff>0)


#delimit ;

heatplot wtp_diff phintBW phintWB if highprob==0, discrete scatter(circle)
colors(blue red, ipolate(20))
p(mlc(black))
keylabels(minmax)
scheme(plotplain)
size(20)
title("Prior<=20%")

;
#delimit cr
graph export "./Graphs/WTP_pattern_low.png", width(1000) height(1000) replace

#delimit ;

heatplot wtp_diff phintBW phintWB if highprob==1, discrete scatter
colors(blue red, ipolate(20))
p(mlc(black))

keylabels(minmax) 
scheme(plotplain)
size(20)
title("Prior>20%")

;
#delimit cr
graph export "./Graphs/WTP_pattern_high.png", width(1000) height(1000) replace


log close


