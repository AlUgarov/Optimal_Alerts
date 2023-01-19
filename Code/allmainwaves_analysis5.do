
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

gen tot_liars=bl_gr+w_gr

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




import delimited using "./Input/coded_strategies.csv", clear


gen mention_hint=0
replace mention_hint=1 if strpos(explanation, "hint")>0
replace mention_hint=1 if strpos(explanation, "says")>0
replace mention_hint=1 if strpos(explanation, "say")>0
replace mention_hint=1 if strpos(explanation, "tell")>0
replace mention_hint=1 if strpos(explanation, "told")>0

gen mention_prop=0

replace mention_prop=1 if strpos(explanation, "balls")

gen mention_hon=0
replace mention_hon=1 if strpos(explanation, "truth")>0
replace mention_hon=1 if strpos(explanation, "honest")>0
replace mention_hon=1 if strpos(explanation, "lying")>0


tab mention_hint
tab mention_prop
tab mention_hon

save "./Temp/coded_strategies.dta", replace
use "./Temp/mainwaves_wide.dta", replace
merge 1:1 participant_id using "./Temp/coded_strategies.dta"
drop _merge



save "./Temp/mainwaves_wide.dta", replace

**BLIND PROTECTION ANALYSIS**
* bp - protection decision (0 - do not protect, 1 - protect)
keep participant_id bp_* bp_time_* sex age stat_educ ncorrect informed_correct wtpq_correct educ
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
collapse (mean) bp submittime (first) sex age stat_educ ncorrect informed_correct wtpq_correct college (sum) totprot=bp (max) switcher backswitcher nbswitches repairable (min) maxspeed=submittime firstswitch=switchround backswitchround, by(participant_id)

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

cluster kmeans ip_w_1 ip_b_1 ip_w_2 ip_b_2 ip_w_3 ip_b_3, k(3)



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

*clean out the subjects with missing informed protection responses (1)
gen miss_ip=missing(ip_w)|missing(ip_b)
bys subject_id: egen incompl_ip=sum(miss_ip)
tab incompl_ip
drop if incompl_ip>0

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
*order tsex sex told old college college tstat_educ stat_educ
*bro


*collapse (mean) sex old stat_educ college final_payoff duration_min (sum) tsex=sex told=old  tcollege=college tstat_educ=stat_educ, by(seq_type)
*order tsex sex told old tcollege college tstat_educ stat_educ
*bro

use "./Output/main_waves.dta", replace

keep subject_id participant_id round time_ip honest ncorrect p phintBW phintWB phintBB rev_response treatm_type tot_liars bl_gr w_gr
collapse (first) treatm_type participant_id time_ip honest ncorrect p phintBW phintWB phintBB rev_response tot_liars bl_gr w_gr, by(subject_id round)
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
label var blackhint "S=Black"
label var honest "Honest treatment"
label var post_prob "Posterior prob."

replace ip_=. if ip_==-99

**Calculate subject-level accuracy of reported beliefs:
gen bel_err=post_prob-be_
label var bel_err "Belief error"

gen absbel_err=abs(bel_err)
sort subject_id
by subject_id: egen tot_bel_err=sum(absbel_err) //total abs error per subject

sum tot_bel_err, detail
local err_med=r(p50)
hist tot_bel_err
gen accur_bel=tot_bel_err<`err_med' //error is less than the median



*Evaluating relative belief accuracy conditional on tasks (prior x signal)
bys round: sum absbel_err
bys subject_id round: egen v1=sum(absbel_err)
bys p phintWB phintBW:egen med_bel_err=median(v1)
gen accur_bel2=abs(bel_err)<med_bel_err

sort  subject_id round hint
order subject_id round phintWB phintBW bel_err absbel_err accur_bel med_bel_err accur_bel2
bro


*Are subjects accurate in the second part of the experiment if they are accurate in the first?
gen laterounds=round>3
bys laterounds subject_id: egen sbel_err=sum(absbel_err)

bys laterounds: sum sbel_err

*How much the signal affects beliefs?
gen be_w=be_
gen be_b=be_
replace be_w=. if hint=="b"
replace be_b=. if hint=="w"
bys subject_id round: egen be_w1=mean(be_w)
bys subject_id round: egen be_b1=mean(be_b)
replace be_w=be_w1
replace be_b=be_b1


gen be_change=be_b-be_w



*Calculate average confidence as reported beliefs that the decision is correct for each signal:
gen post_prob_w=1-post_prob
gen post_prob_b=post_prob
replace post_prob_w=. if hint=="b"
replace post_prob_b=. if hint=="w"
bys subject_id round: egen post_prob_w1=mean(post_prob_w)
bys subject_id round: egen post_prob_b1=mean(post_prob_b)

gen prob_correct=be_
replace prob_correct=1-be_ if hint=="w"
gen prob_B_sign=(p*phintBB+(1-p)*phintBW)
gen prob_W_sign=1-prob_B_sign

gen confid=prob_W_sign*(1-be_w1)+prob_B_sign*be_b1
sum confid
bys plevel phintWB phintBW: sum confid


**Saving belief accuracy:
save "./Temp/base_main_waves.dta", replace
collapse (mean) accur_bel accur_bel2 be_change confid, by(subject_id round)
save "./Temp/bel_accuracy.dta", replace
use "./Temp/base_main_waves.dta", replace



label define accur_bel_l 0 "Inaccurate" 1 "Accur. beliefs"
label values accur_bel accur_bel_l
label values accur_bel2 accur_bel_l

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

merge m:1 participant_id using "./Temp/coded_strategies.dta"
tab strategy_ip
drop _merge

save "./Temp/long_ip_dat.dta", replace
use "./Temp/long_ip_dat.dta", replace
gen xid=1
gen post_probi=round(1000*post_prob)

gen blind_prob=0.05*round

keep if p<0.3

*graph twoway (lpoly bp blind_prob, bwidth(0.1) mlwidth(thick)) (lpoly ip post_prob,  mlwidth(thick)), noscatter title("Protection Responses")

*serrbar ip se post_prob if task_type=="Blind", scale (1.96) title("Blind Protection Response") mlwidth(thick) xtitle("Probability of a black ball") ytitle("Proportion of protection choices")

lpoly ip_ post_prob, bwidth(0.1) ci noscatter title("Informed Protection Response") lineopt(lwidth(1)) mlwidth(thick) xtitle("Posterior probability of a black ball") ytitle("Proportion of protection choices")
graph export "./Graphs/ip_response_lpoly.png", width(1200) height(800) replace

sort subject_id round
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
set more off
****-- INFORMED PROTECTION --****

*expected response by prior prob:
gen highprob=p>0.101
label var highprob "p$=$0.2"
label define highprob_l 0 "p$=$0.1" 1 "p$=$0.2"
label values highprob highprob_l

label var stat_educ "Stat. class"
label define stat_educl 0 "No" 1 "Stat. class"
label values stat_educ stat_educl
tab stat_educ


*Creating splines of prior probability

mkspline post_prob01 0.2 post_prob02 0.4 post_prob03 0.6 post_prob04 0.8 post_prob05 = post_prob

gen high_phintBW=highprob*phintBW
gen high_phintWB=highprob*phintWB
gen black_phintBW=blackhint*phintBW
gen black_phintWB=blackhint*phintWB

gen white_phintBW=(1-blackhint)*phintBW
gen white_phintWB=(1-blackhint)*phintWB

gen FPFN=phintBW*phintWB
gen blackFPFN=phintBW*phintWB
gen whiteFPFN=phintBW*phintWB

label var high_phintBW "FP rate x (p=0.2)"
label var high_phintWB "FN rate x (p=0.2)"
label var black_phintBW "FP rate x (S=Black)"
label var black_phintWB "FN rate x (S=Black)"

label var white_phintBW "FP rate x (S=White)"
label var white_phintWB "FN rate x (S=White)"

label var FPFN "FP rate x FN rate"


xtset subject_id

label var phintBW "FP rate"
label var phintWB "FN rate"

gen ipB=ip_
replace ipB=. if blackhint==0|plevel>201

gen ipW=ip_
replace ipW=. if blackhint==1|plevel>201

bys subject_id: egen varipB=sd(ipB) if plevel<300
bys subject_id: egen varipW=sd(ipW) if plevel<300
sum varipB varipW

sum ip_ if varipB>0&plevel<300&blackhint==1
sum ip_ if varipW>0&plevel<300&blackhint==0
drop varipB varipW ipB ipW

**Baseline IP table:
gen fef=1
eststo clear
eststo: logit ip_ phintBW phintWB blackhint i.plevel if plevel<300, vce(cluster subject_id)
local r2p=e(r2_p)
local llike=e(ll)
test phintBW=phintWB
local p=r(p)
eststo m1: margins, dydx(phintBW phintWB blackhint i.plevel) post
estadd scalar p = `p'
estadd scalar r2p = `r2p'
estadd scalar llike = `llike'

eststo: logit ip_ phintBW phintWB i.plevel if blackhint==0&plevel<300, vce(cluster subject_id)
test phintBW=phintWB
local p=r(p)
local r2p=e(r2_p)
local llike=e(ll)

eststo m2: margins, dydx(phintBW phintWB i.plevel) post
estadd scalar p = `p'
estadd scalar r2p = `r2p'
estadd scalar llike = `llike'

eststo: logit ip_ phintBW phintWB i.plevel if blackhint==1&plevel<300, vce(cluster subject_id)
test phintBW=phintWB
local p=r(p)
local r2p=e(r2_p)
local llike=e(ll)
eststo m3: margins, dydx(phintBW phintWB i.plevel) post
estadd scalar p = `p'
estadd scalar r2p = `r2p'
estadd scalar llike = `llike'

eststo: logit ip_ phintBW phintWB blackhint i.plevel fef i.subject_id if plevel<300, vce(cluster subject_id)
test phintBW=phintWB
local p=r(p)
local r2p=e(r2_p)
local llike=e(ll)
eststo m4: margins, dydx(phintBW phintWB blackhint i.plevel fef) post
estadd scalar p = `p'
estadd scalar r2p = `r2p'
estadd scalar llike = `llike'

eststo: logit ip_ phintBW phintWB i.plevel fef i.subject_id if blackhint==0&plevel<300, vce(cluster subject_id)
test phintBW=phintWB
local p=r(p)
local r2p=e(r2_p)
local llike=e(ll)
eststo m5: margins, dydx(phintBW phintWB i.plevel fef) post
estadd scalar p = `p'
estadd scalar r2p = `r2p'
estadd scalar llike = `llike'

eststo: logit ip_ phintBW phintWB i.plevel fef i.subject_id if blackhint==1&plevel<300, vce(cluster subject_id)
test phintBW=phintWB
local p=r(p)
local r2p=e(r2_p)
local llike=e(ll)
eststo m6: margins, dydx(phintBW phintWB i.plevel fef) post
estadd scalar p = `p'
estadd scalar r2p = `r2p'
estadd scalar llike = `llike'

eststo: logit ip_ phintBW phintWB FPFN i.plevel fef i.subject_id if blackhint==0&plevel<300, vce(cluster subject_id)
test phintBW=phintWB
local p=r(p)
local r2p=e(r2_p)
local llike=e(ll)
eststo m7: margins, dydx(phintBW phintWB FPFN i.plevel fef) post
estadd scalar p = `p'
estadd scalar r2p = `r2p'
estadd scalar llike = `llike'

eststo: logit ip_ phintBW phintWB FPFN i.plevel fef i.subject_id if blackhint==1&plevel<300, vce(cluster subject_id)
test phintBW=phintWB
local p=r(p)
local r2p=e(r2_p)
local llike=e(ll)
eststo m8: margins, dydx(phintBW phintWB FPFN i.plevel fef) post
estadd scalar p = `p'
estadd scalar r2p = `r2p'
estadd scalar llike = `llike'

esttab m1 m2 m3 m4 m5 m6 m7 m8 using "./Tables/table_ip0.tex", b(%9.3g) t(%9.1f) ar2(%9.2f) stats(p N r2p llike, labels("P(FP rate $\neq$ FN rate)" "N" "Pseudo R-squared" "Log-likelihood")) indicate(Subject FE = fef) label addnotes("Errors are clustered by subject, average marginal treatment effects") mtitles("All" "S=White" "S=Black" "All" "S=White" "W=Black" "S=White" "W=Black") title(Informed protection response: logistical regression) star("*" 0.10 "**" 0.05 "***" 0.01) nobaselevels compress nogaps replace


eststo: logit ip_ phintBW phintWB if plevel<300, vce(cluster subject_id)
eststo m2: margins, dydx(phintBW phintWB ) post
eststo: logit ip_ phintBW phintWB blackhint post_prob if plevel<300, vce(cluster subject_id)
eststo m3: margins, dydx(phintBW phintWB blackhint post_prob) post
eststo: logit ip_ phintBW phintWB fef i.subject_id if plevel<300, vce(cluster subject_id)
eststo m4: margins, dydx(phintBW phintWB) post
eststo: logit ip_ phintBW phintWB blackhint fef i.subject_id if plevel<300 , vce(cluster subject_id)
eststo m5: margins, dydx(phintBW phintWB blackhint) post
eststo: logit ip_ phintBW phintWB blackhint post_prob fef i.subject_id if plevel<300, vce(cluster subject_id)
eststo m6: margins, dydx(phintBW phintWB blackhint post_prob) post
esttab m1 m2 m3 m4 m5 m6 using "./Tables/table_ip0x.tex", b(%9.3g) t(%9.1f) ar2(%9.2f) label addnotes("Errors are clustered by subject, average marginal treatment effects") mtitles("" "" "" "FE" "FE" "FE") title(Informed protection response: probit) star("*" 0.10 "**" 0.05 "***" 0.01) nobaselevels compress nogaps replace



eststo clear
eststo: reg ip_ phintBW phintWB  i.plevel if plevel<300, vce(cluster subject_id)
eststo: reg ip_ phintBW phintWB  i.plevel if blackhint==0&plevel<300, vce(cluster subject_id)
eststo: reg ip_ phintBW phintWB i.plevel if blackhint==1&plevel<300, vce(cluster subject_id)
eststo: reg ip_ phintBW phintWB  i.plevel i.subject_id if plevel<300 , vce(cluster subject_id)
eststo: reg ip_ phintBW phintWB  i.plevel i.subject_id if blackhint==0&plevel<300, vce(cluster subject_id)
eststo: reg ip_ phintBW phintWB  i.plevel i.subject_id if blackhint==1&plevel<300, vce(cluster subject_id)
esttab using "./Tables/table_ip0_lin.tex", b(%9.3g) t(%9.1f) ar2(%9.2f) label addnotes("Errors are clustered by subject") indicate(Subject FE = *.subject_id) mtitles("All" "S=White" "S=Black" "All" "S=White" "W=Black") title(Informed protection response: linear regression) star("*" 0.10 "**" 0.05 "***" 0.01) nobaselevels compress nogaps replace




*IP anomalies basic: with flexible control of posteriors:
eststo clear
eststo: logit ip_ post_prob0* phintBW phintWB highprob blackhint, vce(cluster subject_id)
eststo m1: margins, dydx(phintBW phintWB highprob blackhint) post
eststo: logit ip_ post_prob0* phintBW phintWB  highprob blackhint black_phintBW black_phintWB i.subject_id, vce(cluster subject_id)
eststo m2: margins, dydx(phintBW phintWB highprob blackhint black_phintBW black_phintWB) post

eststo: logit ip_ post_prob0* phintBW phintWB  highprob blackhint black_phintBW black_phintWB i.subject_id if plevel<300, vce(cluster subject_id)
eststo m3: margins, dydx(phintBW phintWB highprob blackhint black_phintBW black_phintWB) post

eststo: logit ip_ post_prob0* phintBW phintWB  highprob blackhint black_phintBW black_phintWB i.subject_id if plevel>201, vce(cluster subject_id)
eststo m3: margins, dydx(phintBW phintWB highprob blackhint black_phintBW black_phintWB) post

eststo: logit ip_ post_prob0* phintBW phintWB highprob blackhint high_phintBW high_phintWB i.subject_id if plevel<300, vce(cluster subject_id)
eststo m3: margins, dydx(phintBW phintWB highprob blackhint high_phintBW high_phintWB)  post

eststo: logit ip_ post_prob0* phintBW phintWB highprob blackhint black_phintBW black_phintWB  high_phintBW high_phintWB i.subject_id if plevel<300, vce(cluster subject_id)
eststo m4: margins, dydx(phintBW phintWB  highprob blackhint black_phintBW black_phintWB high_phintBW high_phintWB) post
esttab m1 m2 m3 m4 using "./Tables/table_ip5.tex", b(%9.3g) t(%9.1f) ar2(%9.2f) label addnotes("Reporting average marginal effects, subject FE, errors are clustered by subject." "With flexible controls of posterior probability" ///
 ) mtitles("" "" "" "") title(Informed Protection Response: logit with flexible control for posteriors) star("*" 0.10 "**" 0.05 "***" 0.01) nobaselevels compress nogaps replace

 
*Saving the file to do heterogeneity analysis with respect to IP strategies:
save "./Temp/prep_heter.dta", replace


use "./Temp/prep_heter.dta", replace
*Test average response vs posterior probability by 




**Testing heterogeneity in informed protection response by signal color (black/white) and prior:
set more off
eststo clear
eststo: probit ip_ phintBW phintWB if blackhint==0&plevel==100, vce(cluster subject_id)
eststo: probit ip_ phintBW phintWB if blackhint==0&plevel==200, vce(cluster subject_id)
eststo: probit ip_ phintBW phintWB if blackhint==0&plevel==300, vce(cluster subject_id)
eststo: probit ip_ phintBW phintWB if blackhint==0&plevel==500, vce(cluster subject_id)

eststo: probit ip_ phintBW phintWB  if blackhint==1&plevel==100, vce(cluster subject_id)
eststo: probit ip_ phintBW phintWB  if blackhint==1&plevel==200, vce(cluster subject_id)
eststo: probit ip_ phintBW phintWB  if blackhint==1&plevel==300, vce(cluster subject_id)
eststo: probit ip_ phintBW phintWB  if blackhint==1&plevel==500, vce(cluster subject_id)

esttab using "./Tables/table_ip_het.tex", b(%9.3g) t(%9.1f) ar2(%9.2f) label addnotes("First four for white signal, the rest - black" ///
  "Errors are clustered by subject") mtitles("0.1" "0.2" "0.3" "0.5" "0.1" "0.2" "0.3" "0.5") title(Informed protection by prior) star("*" 0.10 "**" 0.05 "***" 0.01) nobaselevels compress nogaps replace

 
*Generate some fake IP data to see if reacting on the proportion of dishonest gremlins explains our results
gen ip_fake=0
gen false_prob=phintWB+phintBW
replace ip_fake=ip_ if mod(subject_id,3)==0
replace ip_fake=blackhint if mod(subject_id,3)>0
replace ip_fake=(false_prob>0.2) if mod(subject_id,3)==1&hint=="w"
replace ip_fake=(false_prob<0.35) if mod(subject_id,3)==1&hint=="b"
replace ip_fake=(false_prob>0.35) if mod(subject_id,3)==2&hint=="w"

reg ip_fake phintBW phintWB if blackhint==0&plevel<300&mod(subject_id,3)>0

*Is it just an error in beliefs?
*Doing the IP regression and controlling for beliefs:
gen bes=be_^2
mkspline bespline1 0.2 bespline2 0.4 bespline3 0.6 bespline4 0.8 bespline5 = be_

*Studying confidence as a function of false-negative rates: sensitivity seems to slightly decrease with priors
bys plevel blackhint phintWB: sum post_prob
bys plevel: reg post_prob phintWB phintBW if blackhint==1


  
  
**Robustness: flexible control both for beliefs and posteriors:
eststo clear
eststo: logit ip_ post_prob0* bespline* phintBW phintWB if plevel<300, vce(cluster subject_id)
eststo m1: margins, dydx(phintBW phintWB) post
eststo: logit ip_ post_prob0* bespline* phintBW phintWB i.subject_id if plevel<300, vce(cluster subject_id)
eststo m2: margins, dydx(phintBW phintWB) post
eststo: logit ip_ post_prob0* bespline* highprob phintBW phintWB high_phintBW high_phintWB i.subject_id if plevel<300, vce(cluster subject_id)
eststo m3: margins, dydx(highprob phintBW phintWB high_phintBW high_phintWB) post
eststo: logit ip_ post_prob0* bespline* blackhint phintBW phintWB black_phintBW black_phintWB i.subject_id if plevel<300, vce(cluster subject_id)
eststo m4: margins, dydx(blackhint phintBW phintWB black_phintBW black_phintWB) post
eststo: logit ip_ post_prob0* bespline* phintBW phintWB if blackhint==0&plevel<300, vce(cluster subject_id)
eststo m5: margins, dydx(phintBW phintWB) post
eststo: logit ip_ post_prob0* bespline* phintBW phintWB if blackhint==1&plevel<300, vce(cluster subject_id)
eststo m6: margins, dydx(phintBW phintWB) post
esttab m1 m2 m3 m4 m5 m6 using "./Tables/table_ip8_be.tex", b(%9.3g) t(%9.1f) ar2(%9.2f) label addnotes("With flexible controls of posterior probability and beliefs" ///
  "Errors are clustered by subject, average marginal treatment effects") mtitles("" "FE" "" "" "S=White" "S=Black") title(Informed Protection Response: flexible control for posteriors and beliefs) star("*" 0.10 "**" 0.05 "***" 0.01) nobaselevels compress nogaps replace

  
  
**Robustness: flexible control both for beliefs and posteriors:
eststo clear:
eststo: logit ip_ post_prob0* white_phintBW white_phintWB  highprob blackhint black_phintBW black_phintWB i.subject_id  if plevel<300, vce(cluster subject_id)
local r2p=e(r2_p)
local llike=e(ll)
eststo m1: margins, dydx(blackhint white_phintBW white_phintWB black_phintBW black_phintWB highprob) post
estadd scalar r2p = `r2p'
estadd scalar llike = `llike'

eststo: logit ip_ post_prob0* white_phintBW white_phintWB highprob blackhint black_phintBW black_phintWB  high_phintBW high_phintWB i.subject_id if plevel<300, vce(cluster subject_id)
local r2p=e(r2_p)
local llike=e(ll)
eststo m2: margins, dydx(blackhint white_phintBW white_phintWB  black_phintBW black_phintWB highprob high_phintBW high_phintWB) post
estadd scalar r2p = `r2p'
estadd scalar llike = `llike'


eststo: logit ip_ post_prob0* bespline* white_phintBW white_phintWB highprob blackhint black_phintBW black_phintWB i.subject_id if plevel<300, vce(cluster subject_id)
local r2p=e(r2_p)
local llike=e(ll)
eststo m3: margins, dydx(blackhint white_phintBW white_phintWB black_phintBW black_phintWB highprob) post
estadd scalar r2p = `r2p'
estadd scalar llike = `llike'

eststo: logit ip_ post_prob0* bespline* white_phintBW white_phintWB highprob blackhint black_phintBW black_phintWB high_phintBW high_phintWB i.subject_id if plevel<300, vce(cluster subject_id)
local r2p=e(r2_p)
local llike=e(ll)
eststo m4: margins, dydx(blackhint white_phintBW white_phintWB  black_phintBW black_phintWB highprob high_phintBW high_phintWB) post
estadd scalar r2p = `r2p'
estadd scalar llike = `llike'

esttab m1 m2 m3 m4 using "./Tables/table_ip_flex.tex", b(%9.3g) t(%9.1f) label addnotes("With flexible controls of posterior probability and beliefs" ///
  "Subject FE, errors are clustered by subject, average marginal treatment effects") stats(N r2p llike, labels("N" "Pseudo R-squared" "Log-likelihood")) mtitles("Posterior only" "Posterior only" "Both" "Both") title(Informed Protection Response: flexible control for posteriors and beliefs) star("*" 0.10 "**" 0.05 "***" 0.01) nobaselevels compress nogaps replace


  
 **Robustness: LINEAR! flexible control both for beliefs and posteriors:
eststo clear
eststo: reg ip_ post_prob0* white_phintBW white_phintWB  highprob blackhint black_phintBW black_phintWB i.subject_id  if plevel<300, vce(cluster subject_id)
eststo: reg ip_ post_prob0* white_phintBW white_phintWB highprob blackhint black_phintBW black_phintWB  high_phintBW high_phintWB i.subject_id if plevel<300, vce(cluster subject_id)
eststo: reg ip_ post_prob0* bespline* white_phintBW white_phintWB highprob blackhint black_phintBW black_phintWB i.subject_id if plevel<300, vce(cluster subject_id)
eststo: reg ip_ post_prob0* bespline* white_phintBW white_phintWB highprob blackhint black_phintBW black_phintWB high_phintBW high_phintWB i.subject_id if plevel<300, vce(cluster subject_id)
esttab using "./Tables/table_ip_flexlin.tex", b(%9.3g) t(%9.1f)  ar2(%9.2f) label addnotes("With flexible controls of posterior probability and beliefs" ///
  "Subject FE, errors are clustered by subject, average marginal treatment effects") drop(_cons) indicate("Subject FE = *.subject_id" "Posterior=post_prob0*" "Beliefs=bespline*") mtitles("" "" "" "") title(Informed Protection Response: flexible control for posteriors and beliefs, LPM) star("*" 0.10 "**" 0.05 "***" 0.01) nobaselevels compress nogaps replace


  


*Reacting to FP rates is still optimal even when the signal is white, because it lowers the probability of the white signal coming from the white state
* and hence increases the black ball posterior
gen post_prob2=post_prob^2


gen highprobBW=highprob*phintBW
label var highprobBW "FP rate x (p $\geq$ 0.2)"



gen highprobWB=highprob*phintWB
label var highprobWB "FN rate x (p $\geq$ 0.2)"

gen blackhintBW=blackhint*phintBW
label var blackhintBW "FP rate x (S=Black)"
gen blackhintWB=blackhint*phintWB
label var blackhintWB "FN rate x (S=Black)"

gen statBW=stat_educ*phintBW
label var statBW "FP rate x Stat. class"
gen statWB=stat_educ*phintWB
label var statWB "FN rate x Stat. class"

*Final robustness: semiparametric control for posteriors
eststo clear
eststo: semipar ip_ phintBW phintWB if plevel<300, nonpar(post_prob)
eststo: semipar ip_ highprob phintBW highprobBW phintWB highprobWB if plevel<300, nonpar(post_prob) 
eststo: semipar ip_ blackhint phintBW phintWB blackhintBW blackhintWB if plevel<300, nonpar(post_prob)
eststo: semipar ip_ stat_educ phintBW phintWB statBW statWB if plevel<300, nonpar(post_prob)
esttab using "./Tables/table_ip5_semi.tex", b(%9.3g) t(%9.1f) ar2(%9.2f) label mtitles("" "" "" "") title(Informed protection response: semiparametric control for posteriors) star("*" 0.10 "**" 0.05 "***" 0.01) nobaselevels compress nogaps replace





replace bel_err=-bel_err


save "./Temp/prep_beliefs.dta", replace

use "./Temp/prep_beliefs.dta", replace
**Comparison of protection by information 
gen fp_env=phintBW>0
gen fn_env=phintWB>0
bys fp_env fn_env blackhint: sum ip_ if plevel<300

gen prot_cost=5
gen loss=20
gen ip_o=post_prob>(prot_cost/loss)

tempname p1
postfile `p1' false_pos false_neg signal prot ptest posterior mean_opt ptest2 using "./Temp/ip_by_environment.dta", replace

forvalues i=0/1{
  forvalues j=0/1{
  
    forvalues k=0/1{
	  if `k'==0{
	    ttest ip_==0 if plevel<300&fp_env==`i'&fn_env==`j'&blackhint==0, level(95)
		local ptest=r(p_u)
		
	  }
	  else {
	    ttest ip_==1 if plevel<300&fp_env==`i'&fn_env==`j'&blackhint==1, level(95)
	  	local ptest=r(p_l)
	  }
	  local prot=r(mu_1)
	  sum post_prob if plevel<300&fp_env==`i'&fn_env==`j'&blackhint==`k'
	  local posterior=r(mean)
	  ttest ip_==ip_o if plevel<300&fp_env==`i'&fn_env==`j'&blackhint==`k'
	  local mean_opt=r(mu_2)
	  local ptest2=r(p)
	  post `p1' (`i') (`j') (`k') (`prot') (`ptest') (`posterior') (`mean_opt') (`ptest2')
    }
  }
}

postclose `p1'

use "./Temp/ip_by_environment.dta", replace

tostring false_pos false_neg signal, replace
foreach var of varlist false_pos false_neg{
  replace `var'="No" if `var'=="0"
  replace `var'="Yes" if `var'=="1"
}
replace signal="White" if signal=="0"
replace signal="Black" if signal=="1"
bro
format prot ptest posterior mean_opt ptest2 %9.3f
listtex using "./Tables/bigpicture_IP.tex", type rstyle(tabular) head("\begin{table}[H]\centering \footnotesize \caption{Average Protection by Signal Type} \begin{tabular}{cccccccc} \hline \hline" `"\textbf{False-pos.}&\textbf{False-neg.}&\textbf{Signal}&\textbf{\% protect}& \textbf{P(prot$>$0,$<$1)}& \textbf{Posterior} & \textbf{Optimal} & \textbf{P(=optimal)} \\ \hline"') foot("\hline \end{tabular} \end{table}") replace


****-- BELIEF ELICITATION --****

*Beliefs accuracy by signal characteristics:
*bel_err=post_prob-be_

*save "./Temp/prep_beliefs.dta", replace
use "./Temp/prep_beliefs.dta", replace


collapse (mean) ip_ post_prob, by(fp_env fn_env blackhint)
sort fp_env fn_env blackhint

*Anchoring effects for beliefs
*checking if a subject reports the same beliefs for first and second prior:
sort subject_id round blackhint

gen be_p_change=be_~=be_[_n-6]
replace be_p_change=0 if round<4
by subject_id: egen nchanges_p=sum(be_p_change)

gen be_s_change=be_~=be_[_n-2]
replace be_s_change=0 if inlist(round,1,2,4)
by subject_id: egen nchanges_s=sum(be_s_change)
tab nchanges_s nchanges_p



use "./Temp/prep_beliefs.dta", replace
keep if plevel<300
gen fp_env=phintBW>0
gen fn_env=phintWB>0




tempname p1
postfile `p1' false_pos false_neg signal bel_err ptest using "./Temp/bel_by_environment.dta", replace

forvalues i=0/1{
  forvalues j=0/1{
  
    forvalues k=0/1{
	    ttest bel_err==0 if fp_env==`i'&fn_env==`j'&blackhint==`k', level(95)
		local ptest=r(p)
		local bel_err=r(mu_1)
		post `p1' (`i') (`j') (`k') (`bel_err') (`ptest')
	  }
    }
  }

postclose `p1'

use "./Temp/bel_by_environment.dta", replace

tostring false_pos false_neg signal, replace
foreach var of varlist false_pos false_neg{
  replace `var'="No" if `var'=="0"
  replace `var'="Yes" if `var'=="1"
}
replace signal="White" if signal=="0"
replace signal="Black" if signal=="1"
bro
format bel_err ptest %9.3f
listtex using "./Tables/bigpicture_bel.tex", type rstyle(tabular) head("\begin{table}[H]\centering \caption{Average Belief Error by Signal Type} \begin{tabular}{ccccc} \hline \hline" `"\textbf{False-pos.}&\textbf{False-neg.}&\textbf{Signal}&\textbf{Belief error}& \textbf{P($=0$)}\\ \hline"') foot("\hline \end{tabular} \end{table}") replace



*collapse (mean) bel_err, by(fp_env fn_env)
bro
use "./Temp/prep_beliefs.dta", replace

merge m:1 subject_id using "./Temp/ipclasses.dta" //risk aversion, blind prot choices and demographic vars
drop _merge

gen class_alt=2-class
tab class class_alt
label var class_alt "IP strategy class"
label define clsnames 0 "Bayesians" 1 "Simpletons"

label var false_prob "Prop. of lying gremlins"
label values class_alt clsnames


eststo clear
*eststo: reg bel_err phintWB phintBW i.subject_id if plevel<300, vce(cluster subject_id)
*eststo: reg bel_err c.phintWB c.phintBW i.subject_id if plevel<300&blackhint==0, vce(cluster subject_id)
*eststo: reg bel_err c.phintWB c.phintBW i.subject_id if plevel<300&blackhint==1, vce(cluster subject_id)
eststo: reg be_ post_prob blackhint false_prob if plevel<300&class==1, vce(cluster subject_id)
eststo: reg be_ post_prob blackhint false_prob if plevel<300&class==2, vce(cluster subject_id)
esttab using "./Tables/table_be_class.tex", b(%9.3g) se(%9.1f) ar2(%9.2f) label addnotes(Dep. variable: beliefs, errors clustered by subject) title(Belief Elicitation by Class) mtitles("Simpletons" "Cautious Bayesians") star("*" 0.10 "**" 0.05 "***" 0.01) nobaselevels compress nogaps replace




eststo clear
eststo: reg bel_err phintWB phintBW i.subject_id if plevel<300, vce(cluster subject_id)
eststo: reg bel_err phintWB phintBW i.subject_id if plevel<300&blackhint==0, vce(cluster subject_id)
eststo: reg bel_err phintWB phintBW i.subject_id if plevel<300&blackhint==1, vce(cluster subject_id)
esttab using "./Tables/table_be_err.tex", b(%9.3g) se(%9.1f) ar2(%9.2f) label addnotes(Dep. variable: reported belief - posterior probability) title(Belief Elicitation: When Mistakes Happen) mtitles("All" "S=White" "S=Black") star("*" 0.10 "**" 0.05 "***" 0.01) indicate(Subject FE = *.subject_id) nobaselevels compress nogaps replace


eststo clear
eststo: reg bel_err i.class_alt##c.phintWB i.class_alt##c.phintBW i.subject_id if plevel<300, vce(cluster subject_id)
eststo: reg bel_err i.class_alt##c.phintWB i.class_alt##c.phintBW i.subject_id if plevel<300&blackhint==0, vce(cluster subject_id)
eststo: reg bel_err i.class_alt##c.phintWB i.class_alt##c.phintBW i.subject_id if plevel<300&blackhint==1, vce(cluster subject_id)
esttab using "./Tables/table_be_errx.tex", b(%9.3g) se(%9.1f) ar2(%9.2f) label addnotes(Dep. variable: reported belief - posterior probability) title(Belief Elicitation: When Mistakes Happen) mtitles("All" "S=White" "S=Black") star("*" 0.10 "**" 0.05 "***" 0.01) indicate(Subject FE = *.subject_id) nobaselevels compress nogaps replace






reg bel_err i.class_alt##c.phintWB i.class_alt##c.phintBW i.subject_id if plevel<300&blackhint==0, vce(cluster subject_id)


*Testing the accuracy of reported beliefs***
eststo clear
eststo: reg bel_err phintWB phintBW if plevel<300, vce(cluster subject_id)
eststo: reg bel_err i.plevel phintWB phintBW if plevel<300, vce(cluster subject_id)
eststo: reg bel_err i.goodquiz##c.phintWB i.goodquiz##c.phintBW if plevel<300, vce(cluster subject_id)
eststo: reg bel_err i.goodquiz##i.plevel i.goodquiz##c.phintWB i.goodquiz##c.phintBW if plevel<300, vce(cluster subject_id)
eststo: reg bel_err i.stat_educ##c.phintWB i.stat_educ##c.phintBW if plevel<300, vce(cluster subject_id)
eststo: reg bel_err i.stat_educ##i.plevel i.stat_educ##c.phintWB i.stat_educ##c.phintBW if plevel<300, vce(cluster subject_id)
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
eststo: reg lt_bel lt_prior signal if plevel<300, noconstant vce(cluster subject_id)
eststo: xtreg lt_bel lt_prior signal if plevel<300, fe vce(cluster subject_id)
eststo: reg lt_bel lt_prior signal i.goodquiz#c.lt_prior i.goodquiz#c.signal if plevel<300, noconstant vce(cluster subject_id)
eststo: xtreg lt_bel lt_prior signal i.goodquiz#c.lt_prior i.goodquiz#c.signal if plevel<300, fe vce(cluster subject_id)
eststo: reg lt_bel lt_prior signal i.stat_educ#c.lt_prior i.stat_educ#c.signal if plevel<300, noconstant vce(cluster subject_id)
eststo: xtreg lt_bel lt_prior signal i.stat_educ#c.lt_prior i.stat_educ#c.signal if plevel<300, fe vce(cluster subject_id)
esttab using "./Tables/table_be3.tex", b(%9.3g) t(%9.1f) ar2(%9.2f) label drop(_cons) mtitles("OLS" "FE" "OLS" "FE" "OLS" "FE") star("*" 0.10 "**" 0.05 "***" 0.01) note("Decomposition works only for imperfect signals") nobaselevels compress nogaps replace

gen psignalblack=p*phintBB+(1-p)*phintBW
gen psignalwhite=1-psignalblack
gen pifsignblack=p*phintBB/psignalblack


**********************************************
****-- Main Regressions: WTP FOR INFORMATION --****
**********************************************
use "./Output/main_waves.dta", replace
xtset subject_id round

merge m:1 subject_id round using "./Temp/bel_accuracy.dta" //average belief accuracy by subject
drop _merge


merge m:1 participant_id p using "./Temp/bp_val.dta" //risk aversion, blind prot choices and demographic vars
drop if _merge==2
drop _merge

merge m:1 subject_id using "./Temp/ipclasses.dta" //risk aversion, blind prot choices and demographic vars
drop _merge

drop if pilot==1 //dropping the pilot

gen fp_env=phintBW>0
gen fn_env=phintWB>0


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
label value accur_bel2 accur_bel_l


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


scatter wtp value  if p<0.3, title("Theoretical vs actual WTP") xtitle("Theoretical WTP") ytitle("Actual WTP") jitter(2) 
graph export "./Graphs/WTP_value_scatter.png", width(1200) height(800) replace

heatplot wtp value  if p<0.3, normalize backfill colors(HCL blues, reverse) title("Theoretical vs actual WTP") xtitle("Theoretical WTP") ytitle("Actual WTP") levels(10) bins(10) keylabels(minmax)
graph export "./Graphs/WTP_value_heat.png", width(1200) height(800) replace

heatplot wtp value_ra  if p<0.3, normalize backfill colors(HCL blues, reverse) title("Theoretical with risk aversion vs actual WTP") xtitle("Theoretical WTP") ytitle("Actual WTP") levels(10) bins(10) keylabels(minmax)
graph export "./Graphs/WTP_value_heat_ra.png", width(1200) height(800) replace

*Note: correlation wtp value - 0.27; correlation wtp value_ra - 0.24
heatplot wtp value  if p<0.3, colors(s2) title("Theoretical vs actual WTP") xtitle("Theoretical WTP") ytitle("Actual WTP") discrete levels(1) scatter sizeprop keylabels(none)
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

replace accur_bel2=1-round(accur_bel2)
label value accur_bel2 accur_bel1l

*Self-reported strategies
*gen strategy_short=strategy_ip
*replace strategy_short=3 if strategy_short>3
label var strategy_short "Strategy"
label define strategy_short_l 1 "Rational" 2 "Seek honest" 3 "Other"
label values strategy_short strategy_short_l


label var be_change "Belief change"
label var confid "Certainty"



bys fp_env fn_env: ttest wtp_diff == 0 if plevel<300, level(95)



**Baseline WTP regs: tobit
eststo clear
eststo: tobit wtp false_neg false_pos if plevel<300, ll(0) ul(5)
eststo: tobit wtp false_neg false_pos cost_bp if plevel<300, ll(0) ul(5)
eststo: tobit wtp i.risk_pref##c.false_pos i.risk_pref##c.false_neg if plevel<300, ll(0) ul(5)
eststo: tobit wtp i.accur_bel2##c.false_pos i.accur_bel2##c.false_neg if plevel<300, ll(0) ul(5)
eststo: tobit wtp i.plevel##c.false_pos i.plevel##c.false_neg if plevel<300, ll(0) ul(5)
esttab using "./Tables/table_wtp_01tob.tex", b(%9.3g) se(%9.1f) ar2(%9.2f) label title(WTP for Information (tobit)) mtitles("" "" "" "" "") star("*" 0.10 "**" 0.05 "***" 0.01) nobaselevels compress nogaps replace


**Baseline WTP regs: tobit

gen class_alt=2-class

tab class class_alt
label var class_alt "IP strategy class"
label define clsnames 0 "Bayesians" 1 "Simpletons"
label values class_alt clsnames


eststo clear
eststo: tobit wtp false_neg false_pos if plevel<300, ll(0) ul(5)
eststo: tobit wtp i.class_alt##c.false_neg i.class_alt##c.false_pos if plevel<300, ll(0) ul(5)
eststo: tobit wtp false_neg false_pos, ll(0) ul(5)
eststo: tobit wtp i.class_alt##c.false_neg i.class_alt##c.false_pos, ll(0) ul(5)
esttab using "./Tables/table_wtp_classes.tex", b(%9.3g) se(%9.1f) ar2(%9.2f) label title(WTP for Information: heterogeneity by IP class) mtitles("p$<$0.3" "p$<$0.3" "All" "All") star("*" 0.10 "**" 0.05 "***" 0.01) nobaselevels compress nogaps replace





**Baseline WTP regs 2: tobit
eststo clear
eststo: tobit wtp false_neg false_pos if plevel<300, ll(0) ul(5)
eststo: tobit wtp false_neg false_pos if plevel==100, ll(0) ul(5)
eststo: tobit wtp false_neg false_pos if plevel==200, ll(0) ul(5)
eststo: tobit wtp false_neg false_pos cost_bp if plevel<300, ll(0) ul(5)
eststo: tobit wtp false_neg false_pos cost_bp be_change if plevel<300, ll(0) ul(5)
eststo: tobit wtp false_neg false_pos cost_bp confid if plevel<300, ll(0) ul(5)
esttab using "./Tables/table_wtp_02tob.tex", b(%9.3g) se(%9.1f) ar2(%9.2f) label title(WTP for Information (tobit)) mtitles("All" "p=0.1" "p=0.2" "All" "All" "All") star("*" 0.10 "**" 0.05 "***" 0.01) nobaselevels compress nogaps replace


**Strategy effects
eststo clear
eststo: reg wtp_diff i.strategy_short false_neg false_pos if plevel<300
eststo: reg wtp_diff i.strategy_short##c.false_neg i.strategy_short##c.false_pos if plevel<300
eststo: reg wtp false_neg false_pos if plevel==100
eststo: reg wtp_diff i.strategy_short##c.false_neg i.strategy_short##c.false_pos if plevel==100
eststo: reg wtp false_neg false_pos if plevel==200
eststo: reg wtp_diff i.strategy_short##c.false_neg i.strategy_short##c.false_pos if plevel==200
esttab using "./Tables/table_wtp_strat.tex", b(%9.3g) se(%9.1f) ar2(%9.2f) label title(WTP minus Value of Information, connection to self-reported protection strategy) mtitles("All" "p=0.1" "p=0.2" "All" "All" "All") star("*" 0.10 "**" 0.05 "***" 0.01) nobaselevels compress nogaps replace


eststo clear
eststo: reg wtp_diff i.strategy_short  phintBW phintWB if plevel<300
eststo: reg wtp_diff i.strategy_short##c.phintBW i.strategy_short##c.phintWB if plevel<300
eststo: reg wtp phintBW phintWB if plevel==100
eststo: reg wtp_diff i.strategy_short##c.phintBW i.strategy_short##c.phintWB if plevel==100
eststo: reg wtp phintBW phintWB if plevel==200
eststo: reg wtp_diff i.strategy_short##c.phintBW i.strategy_short##c.phintWB if plevel==200



**Baseline WTP difference: risk-aversion and belief accuracy:
eststo clear
eststo: reg wtp_diff false_pos false_neg if plevel<300, vce(cluster subject_id)
eststo: reghdfe wtp_diff false_pos false_neg if plevel<300, abs(subject_id) vce(cluster subject_id)
eststo: reghdfe wtp_diff i.risk_pref##c.false_pos i.risk_pref##c.false_neg if plevel<300, abs(subject_id) vce(cluster subject_id)
eststo: reghdfe wtp_diff i.accur_bel2##c.false_pos i.accur_bel2##c.false_neg if plevel<300, abs(subject_id) vce(cluster subject_id)
eststo: reghdfe wtp_diff i.plevel##c.false_pos i.plevel##c.false_neg if plevel<300, abs(subject_id) vce(cluster subject_id)
*eststo: reghdfe wtp_diff i.accur_bel2#i.risk_pref##c.false_pos i.accur_bel2#i.risk_pref##c.false_neg if plevel<300, abs(subject_id) vce(cluster subject_id)
esttab using "./Tables/table_wtpdiff_01rfull_ag.tex", b(%9.3g) se(%9.1f) ar2(%9.2f) label title(WTP minus Value of Information (OLS)) mtitles("" "" "" "" "") star("*" 0.10 "**" 0.05 "***" 0.01) nobaselevels compress nogaps replace
esttab using "./Tables/table_wtpdiff_01r_ag.tex", b(%9.3g) se(%9.1f) ar2(%9.2f) label title(WTP minus Value of Information (OLS)) mtitles("" "" "" "" "") star("*" 0.10 "**" 0.05 "***" 0.01) noomit nobaselevels compress nogaps replace





*Demographic variables:
eststo clear
eststo: reg wtp_diff false_pos false_neg  if plevel<300, vce(cluster subject_id)
eststo: reg wtp_diff i.sex##c.false_pos i.sex##c.false_neg  if plevel<300, vce(cluster subject_id)
eststo: reg wtp_diff i.sex##i.plevel i.sex##c.false_pos i.sex##c.false_neg if plevel<300, vce(cluster subject_id)
eststo: reg wtp_diff i.stat_educ##c.false_pos i.stat_educ##c.false_neg if plevel<300, vce(cluster subject_id)
eststo: reg wtp_diff i.stat_educ##i.plevel i.stat_educ##c.false_pos i.stat_educ##c.false_neg if plevel<300, vce(cluster subject_id)
eststo: reg wtp_diff i.old##c.false_pos i.old##c.false_neg if plevel<300, vce(cluster subject_id)
eststo: reg wtp_diff i.old##i.plevel i.old##c.false_pos i.old##c.false_neg if plevel<300, vce(cluster subject_id)
eststo: reg wtp_diff i.goodquiz##c.false_pos i.goodquiz##c.false_neg if plevel<300, vce(cluster subject_id)
eststo: reg wtp_diff i.goodquiz##i.plevel i.goodquiz##c.false_pos i.goodquiz##c.false_neg if plevel<300, vce(cluster subject_id)
esttab using "./Tables/table_wtpdiff_02.tex", b(%9.3g) se(%9.1f) ar2(%9.2f) label indicate(Prior dummies=*.plevel) title(WTP minus Value of Information: demographic determinants) mtitles("" "" "" "" "" "" "" "" "") star("*" 0.10 "**" 0.05 "***" 0.01) nobaselevels compress nogaps replace




*By prior prob:
eststo clear
eststo: reghdfe wtp_diff false_pos false_neg if plevel==100, abs(subject_id) vce(cluster subject_id)
eststo: reghdfe wtp_diff false_pos false_neg if plevel==200, abs(subject_id) vce(cluster subject_id)
eststo: reghdfe wtp_diff false_pos false_neg if plevel==300, abs(subject_id) vce(cluster subject_id)
eststo: reghdfe wtp_diff false_pos false_neg if plevel==500, abs(subject_id) vce(cluster subject_id)
esttab using "./Tables/table_wtpdiff_03.tex", b(%9.3g) se(%9.1f) ar2(%9.2f) label title(WTP - Value of Information, by prior) addnotes("Subject fixed effects are included.") mtitles("0.1" "0.2" "0.3" "0.5") star("*" 0.10 "**" 0.05 "***" 0.01) nobaselevels compress nogaps replace



*By prior prob:
eststo clear
eststo: tobit wtp phintBW phintWB if plevel==100, ll(0) ul(5)
eststo: tobit wtp phintBW phintWB if plevel==200, ll(0) ul(5)
eststo: tobit wtp phintBW phintWB if plevel==300, ll(0) ul(5)
eststo: tobit wtp phintBW phintWB if plevel==500, ll(0) ul(5)
esttab using "./Tables/table_wtpdiff_04tob.tex", b(%9.3g) se(%9.1f) ar2(%9.2f) label title(WTP for Information, by prior (tobit)) mtitles("0.1" "0.2" "0.3" "0.5") star("*" 0.10 "**" 0.05 "***" 0.01) nobaselevels compress nogaps replace


save  "./Temp/wtp_discrepancy0.dta", replace

use "./Temp/wtp_discrepancy0.dta", replace

**Calculate average WTP discrepancy by signal type
tempname p1
postfile `p1' false_pos false_neg wtp_diff ptest using "./Temp/wtp_by_environment.dta", replace

forvalues i=0/1{
  forvalues j=0/1{
	    ttest wtp_diff==0 if fp_env==`i'&fn_env==`j'&plevel<300, level(95)
		local ptest=r(p)
		local wtp_diff=r(mu_1)
		post `p1' (`i') (`j') (`wtp_diff') (`ptest')
    }
  }

postclose `p1'

use "./Temp/wtp_by_environment.dta", replace

tostring false_pos false_neg, replace
foreach var of varlist false_pos false_neg{
  replace `var'="No" if `var'=="0"
  replace `var'="Yes" if `var'=="1"
}
bro
format wtp_diff ptest %9.3f
listtex using "./Tables/bigpicture_wtp.tex", type rstyle(tabular) head("\begin{table}[H]\centering \caption{Average WTP discrepancy (WTP-Value) by Signal Type} \begin{tabular}{cccc} \hline \hline" `"\textbf{False-positive}&\textbf{False-negative}&\textbf{Mean WTP discrepancy}& \textbf{P($=0$)}\\ \hline"') foot("\hline \end{tabular} \end{table}") replace





*How many switch choices after changing prior probabilities:
gen wtp_p_change=wtp~=wtp[_n-3]
replace wtp_p_change=0 if round<4
by subject_id: egen nchanges=sum(wtp_p_change)

gen wtp_s_change=wtp~=wtp[_n-1]
replace wtp_s_change=0 if inlist(round,1,2,4)
by subject_id: egen nchanges_s=sum(wtp_s_change)
tab nchanges_s

**Removing subjects which do not update their WTP after increasing priors:
eststo clear
eststo: reghdfe wtp_diff phintBW phintWB if nchanges>0, abs(subject_id) vce(cluster subject_id)
eststo: reghdfe wtp_diff phintBW phintWB if plevel==100&nchanges>0, abs(subject_id) vce(cluster subject_id)
eststo: reghdfe wtp_diff phintBW phintWB if plevel==200&nchanges>0, abs(subject_id) vce(cluster subject_id)
eststo: reghdfe wtp_diff phintBW phintWB if plevel==300&nchanges>0, abs(subject_id) vce(cluster subject_id)
eststo: reghdfe wtp_diff phintBW phintWB if plevel==500&nchanges>0, abs(subject_id) vce(cluster subject_id)
esttab using "./Tables/table_wtpdiff_06s.tex", b(%9.3g) se(%9.1f) ar2(%9.2f) label title(WTP - Value of Information, by prior) mtitles("All" "0.1" "0.2" "0.3" "0.5")  addnotes("Only subjects who change their decisions across priors") star("*" 0.10 "**" 0.05 "***" 0.01) nobaselevels compress nogaps replace


*Do belief changes matter? be_change is the difference in reported beliefs between black and white signals
bys plevel: reg wtp_diff phintBW phintWB be_change, vce(cluster subject_id)

bys plevel: reg wtp_diff phintBW phintWB confid, vce(cluster subject_id)

*By prior prob with high/low prob: testing for the presentation order effect
*No variation in the presentation order. Comparing 0.1 vs 0.2 and 0.3 vs 0.5
gen highsession=seq>3
gen firstorder=round<4

label var highsession "Higher probabilities session"
label define highsessionl 0 "Low probabilities" 1 "Starts with p=0.2"
label value highsession highsessionl

label var firstorder "First prior in the sequence"
label define firstorderl 0 "Second prior" 1 "First prior"
label value firstorder firstorderl

tab highsession firstorder


eststo clear
eststo: reg wtp_diff phintBW phintWB if p<0.3, vce(cluster subject_id)
eststo: reg wtp_diff phintBW phintWB if p>0.25, vce(cluster subject_id)
eststo: reg wtp_diff i.highsession##c.phintBW i.highsession##c.phintWB if p<0.3, vce(cluster subject_id)
*eststo: reg wtp_diff phintBW phintWB if highsession==0, vce(cluster subject_id)
*eststo: reg wtp_diff phintBW phintWB if highsession==1, vce(cluster subject_id)
eststo: reg wtp_diff i.highsession##c.phintBW i.highsession##c.phintWB, vce(cluster subject_id)
eststo: reg wtp_diff i.firstorder##c.phintBW i.firstorder##c.phintWB, vce(cluster subject_id)
eststo: reg wtp_diff i.firstorder##c.phintBW i.firstorder##c.phintWB i.highsession#c.phintBW i.highsession#c.phintWB, vce(cluster subject_id)
esttab using "./Tables/table_wtpdiff_order.tex", b(%9.3g) se(%9.1f) ar2(%9.2f) label title(WTP - Value of Information, by prior with order effects) mtitles( "p=0.1,0.2" "p=0.3,0.5" "p=0.1,0.2" "" "" "") star("*" 0.10 "**" 0.05 "***" 0.01) nobaselevels compress nogaps replace






*Testing heterogeneity
reg wtp_diff i.plevel##c.false_pos i.plevel##c.false_neg, vce(cluster subject_id)
contrast plevel plevel#c.false_pos plevel#c.false_neg, overall


*By prior prob 2:
eststo clear
eststo: reghdfe  wtp_diff phintBW phintWB, abs(subject_id) vce(cluster subject_id)
eststo: reghdfe  wtp_diff phintBW phintWB if plevel==100, abs(subject_id) vce(cluster subject_id)
eststo: reghdfe  wtp_diff phintBW phintWB if plevel==200, abs(subject_id) vce(cluster subject_id)
eststo: reghdfe  wtp_diff phintBW phintWB if plevel==300, abs(subject_id) vce(cluster subject_id)
eststo: reghdfe  wtp_diff phintBW phintWB if plevel==500, abs(subject_id) vce(cluster subject_id)
esttab using "./Tables/table_wtpdiff_06.tex", b(%9.3g) se(%9.1f) ar2(%9.2f) label title(WTP - Value of Information, by prior) mtitles("All" "0.1" "0.2" "0.3" "0.5") star("*" 0.10 "**" 0.05 "***" 0.01) nobaselevels compress nogaps replace



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
eststo: tobit wtp i.accur_bel1##c.highprob i.accur_bel1##c.phintBW i.accur_bel1##i.highprob#c.phintBW i.accur_bel1##c.phintWB i.accur_bel1##i.highprob#c.phintWB, ll(0) ul(5)
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


save "./Temp/wtpdat2.dta", replace
use "./Temp/wtpdat2.dta", replace
collapse (mean) wtp value (count) nobs=value, by(plevel phintBW phintWB)

gen phintBWs=round(10*phintBW)
gen phintWBs=round(10*phintWB)
gen plevels=round(0.01*plevel)
tostring phintBWs phintWBs plevels, force replace
gen treatm_type=plevels+" "+phintBWs+" "+phintWBs
tab treatm_type
twoway (scatter wtp value, mlabel(treatm_type) msize(large)) ( function y = x, ra(0 3.5)) if plevel<300, title("Value vs average WTP by treatment") note("Treatment labels: first number - probability, then FN, FP rates. Low priors only.") ylabel(0(0.5)3.5) xlabel(0(0.5)3.5) legend(off)

graph export "./Graphs/WTP_value_average.png", width(1200) height(800) replace


use "./Temp/wtpdat2.dta", replace

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


