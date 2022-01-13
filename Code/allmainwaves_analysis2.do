
**ANALYZE THE FIRST WAVE (planned 100 participants)***
set more off
clear all

**Requires: estout, moremata

*!!put your working folder below:

cd C:\Tornado_warnings\Experiment\Alerts_Experiment



log using "./Temp/pilot_analysis.log", replace
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

*sex variable:
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


***SANITY CHECKS*******************

**??**

*encode participant_id, gen(subject_id) //as participant_id is initially a string



**BLIND PROTECTION ANALYSIS**
* bp - protection decision (0 - do not protect, 1 - protect)

keep participant_id bp_* bp_time_* sex age stat_educ ncorrect
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


*Collapsing to have one obs per participant (study participant's characteristics):
collapse (mean) bp submittime (first) sex age stat_educ ncorrect (sum) totprot=bp (max) switcher backswitcher nbswitches repairable (min) maxspeed=submittime firstswitch=switchround backswitchround, by(participant_id)

tab firstswitch
replace firstswitch=7-totprot if repairable==1
tab firstswitch

gen backswitcher0=backswitcher
replace backswitcher=0 if repairable==1

tab switcher
tab backswitcher
tab firstswitch //the distribution of the first switching round
tab backswitchround

//scatter submittime backswitcher

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

append using "./Temp/bp_val_pilot.dta"
replace pilot=1 if missing(pilot)

keep participant_id p bp bp_val pilot

save "./Temp/bp_val.dta", replace













**CHANGE TO THE PANEL STRUCTURE FOR THE OTHER TASKS********************
**Panel: participant_id round 
**because all the tasks except the blind protection have the same N of rounds (6) 
use "./Temp/mainwaves_wide.dta", replace

*add the pilot's data:
append using "./Temp/pilot_wide.dta"
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


//Calculate exp costs under the optimal strategy:
gen ip_val_o=-(p*(phintWB*(1-ip_w_o)+phintBB*(1-ip_b_o))*loss+p*(phintWB*ip_w_o+phintBB*ip_b_o)*protectioncost+(1-p)*(phintWW*ip_w_o+phintBW*ip_b_o)*protectioncost)


label var ip_val "Exp. costs"
label var ip_val_o "Optimal exp. costs"


gen ip_val_diff=ip_val-ip_val_o //discrepancy between actual and optimal expected costs of informed protection

list subject_id p ip_w ip_b ip_val ip_val_o ip_val_diff if abs(ip_val_diff)>3
//large discrepancies emerge when subjects take contrarian actions in the informed protection (e.g. protect when white and do not protect when black)


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





*Saving the cleaned dataset with the panel structure
save "./Output/main_waves.dta", replace

keep subject_id participant_id round time_ip honest ncorrect p phintBW phintWB phintBB rev_response
collapse (first) participant_id time_ip honest ncorrect p phintBW phintWB phintBB rev_response, by(subject_id round)
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


**********************************************
****----MAIN REGRESSIONS---***********
**********************************************
xtset subject_id question


****-- INFORMED PROTECTION --****

**response to posterior probabilities of black***
eststo clear
eststo: probit ip_ post_prob, vce(robust)
eststo: probit ip_ post_prob p blackhint, vce(robust)
eststo: probit ip_ post_prob if ncorrect>6, vce(robust)
eststo: probit ip_ post_prob p blackhint if ncorrect>6, vce(robust)
esttab using "./Tables/table_ip1.tex", b(%9.3g) t(%9.1f) aic(%9.2f) label title(Informed Protection) mtitles("All" "All" "Good quiz" "Good quiz") star("*" 0.10 "**" 0.05 "***" 0.01) compress nogaps replace



**response to elicited (and posterior) probabilities of black***
eststo clear
eststo: probit ip_ be_, vce(robust)
eststo: probit ip_ be_ bel_err, vce(robust)
eststo: probit ip_ be_ bel_err if ncorrect>6, vce(robust)
esttab using "./Tables/table_ip2.tex", b(%9.3g) t(%9.1f) aic(%9.2f) label title(Informed Protection: Response to Reported Beliefs) mtitles("All" "All" "Good quiz") star("*" 0.10 "**" 0.05 "***" 0.01) compress nogaps replace



eststo clear
eststo: probit ip_ be_, vce(robust)
eststo: probit ip_ be_ bel_err, vce(robust)
eststo: probit ip_ i.goodquiz##c.be_ i.goodquiz##c.bel_err, vce(robust)
eststo: probit ip_ i.stat_educ##c.be_ i.stat_educ##c.bel_err, vce(robust)
esttab using "./Tables/table_ip3.tex", b(%9.3g) t(%9.1f) aic(%9.2f) label title(Informed Protection: Response to Reported Beliefs) mtitles("" "" "" "" "" "") star("*" 0.10 "**" 0.05 "***" 0.01) nobaselevels compress nogaps replace



**testing for timing and correlation
xtreg ip_ be_, fe vce(robust)
xtreg ip_ be_ if time_ip<20, fe vce(robust)
reg ip_ post_prob, vce(robust)
reg ip_ post_prob if time_ip<20, vce(robust)



****-- BELIEF ELICITATION --****
**prepare variables for belief updating responsiveness analysis (Mobius et al, 2011)
**our sample is small because sme vars req to calculate log(0)
gen lt_bel=log(be_/(1-be_))
gen lt_prior=log(p/(1-p))
gen lambda_B=log(phintBB/phintBW)
gen lambda_W=log((1-phintBB)/(1-phintBW))
gen lt_posterior=log(post_prob/(1-post_prob)) //just for verification
gen signalB=blackhint*lambda_B
gen signalW=(1-blackhint)*lambda_W
sum lt_bel lt_posterior lt_prior signalB signalW

reg lt_posterior lt_prior signalB signalW //all the coeffs should be one, zero const



*Testing the accuracy of reported beliefs***
eststo clear
eststo: reg be_ post_prob, vce(robust)
eststo: reg be_ post_prob if ncorrect>6, vce(robust)
eststo: reg be_ post_prob if (honest==0), vce(robust)
eststo: reg be_ post_prob if (accur_bel==1), vce(robust)
esttab using "./Tables/table_be1.tex", b(%9.3g) t(%9.1f) ar2(%9.2f) label title(Belief Elicitation: Belief vs Posterior) mtitles("All" "Good quiz" "Dishonest greml" "Acc. beliefs") star("*" 0.10 "**" 0.05 "***" 0.01) compress nogaps replace



eststo clear
eststo: reg be_ post_prob, vce(robust)
test post_prob=1
estadd scalar CoefVarName = r(F)
eststo: reg be_ post_prob if ncorrect>6, vce(robust)
test post_prob=1
estadd scalar CoefVarName = r(F)
eststo: reg be_ post_prob if (honest==0), vce(robust)
test post_prob=1
estadd scalar CoefVarName = r(F)
esttab using "./Tables/table_be1.tex", b(%9.3g) t(%9.1f) ar2(%9.2f) stats() label title(Belief Elicitation: Belief vs Posterior) mtitles("All" "Not_honest" "Good quiz") star("*" 0.10 "**" 0.05 "***" 0.01) compress nogaps replace




*Decomposition of belief updating (coeffs should be all ones), small sample***
eststo clear
eststo: reg lt_bel lt_prior signalB signalW, noconstant vce(robust)
eststo: xtreg lt_bel lt_prior signalB signalW, fe vce(robust)
eststo: xtreg lt_bel lt_prior signalB signalW if ncorrect>6, fe vce(robust)
esttab using "./Tables/table_be3.tex", b(%9.3g) t(%9.1f) ar2(%9.2f) label title(Belief Elicitation: Decomposition) mtitles("OLS" "FE" "Good quiz, FE") star("*" 0.10 "**" 0.05 "***" 0.01) compress nogaps replace


****MEASURING THE INFORMED PROTECTION RESPONSE TO BELIEFS VS POSTERIORS AND GRAPHING IT:
**Calculate the N of reversions in strategy with respect to 1) posterior probab and 2) beliefs:
sort subject_id ip_
by subject_id: egen rpostprob=rank(post_prob)
gen rpostprob0=(1-ip_)*rpostprob
gen rpostprob1=ip_*rpostprob

by subject_id: egen rbe=rank(be_)
gen rbe0=(1-ip_)*rbe
gen rbe1=ip_*rbe

gen bel_errsq=bel_err^2

collapse (first) tot_bel_err accur_bel (sum) bel_errsq sumrank_prob0=rpostprob0 sumrank_prob1=rpostprob1 sumrank_bel0=rbe0 sumrank_bel1=rbe1 n1=ip_, by(subject_id)
label var accur_bel "Accur. beliefs"
label values accur_bel accur_bel_l

gen correction=n1*(n1+1)/2

sum sumrank_prob0 sumrank_prob1
gen U_postprob=sumrank_prob1-correction
gen U_bel=sumrank_bel1-correction

sum U_postprob U_bel

gen bel_precision=1/sqrt(bel_errsq)

gen special=U_bel>20&bel_precision<1
replace special=special+1

gen bel_precisionw=sqrt(bel_precision) //a bit of scaling to make graphs prettier

scatter U_postprob U_bel [w=bel_precisionw], title("Determinants of informed protection response (by subject)") xtitle("Mann-Whitney U (beliefs)") ytitle("Mann-Whitney U (posteriors)") note("Markers' area is inv.proportional to beliefs accuracy (sqrt(sum(error^2))")
graph export "./Graphs/clustering.png", width(1200) height(800) replace


scatter U_bel bel_precision, title("Determinants of informed protection response (by subject)") xtitle("Beliefs' accuracy") ytitle("Mann-Whitney U (beliefs)") 
graph export "./Graphs/clustering2.png", width(1200) height(800) replace

//there are respondents with inaccurate beliefs who still follow them, but they also follow posteriors:

//in short, people with inaccurate beliefs either behave inconsistently with everything or behave consistently with posteriors:
sum tot_bel_err, detail
list U_postprob U_bel bel_precision tot_bel_err if U_bel>20&bel_precision<1


save "./Temp/bel_accuracy.dta", replace







****-- WTP FOR INFORMATION --****
use "./Output/main_waves.dta", replace
xtset subject_id round

merge m:1 subject_id using "./Temp/bel_accuracy.dta" //average belief accuracy by subject
drop _merge


merge m:1 participant_id p using "./Temp/bp_val.dta" //risk aversion, blind prot choices and demographic vars
drop if _merge==2
drop _merge

*calculate theoretical value of signal for a risk-neutral subject:
gen cost_bp=min(p*loss, protectioncost) //expected blind protection cost
gen false_pos=(1-p)*phintBW*protectioncost //expected protection costs for false positives
gen true_pos=p*phintBB*protectioncost //expected protection costs for true positives



gen false_neg=p*phintWB*loss //expected false positive costs
gen cost_ip=false_neg+false_pos+true_pos //total informed protection costs
gen pos_cost=false_pos+true_pos //total protection costs in response to positive signals
gen value=max(0, cost_bp-cost_ip) //theoretical value for a risk-neutral subject

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
  X = st_data(.,("p","phintWW","phintBB","theta"))
  Z=J(rows(X),1,0)
  
  Z1=J(rows(X),1,0)
  Z2=J(rows(X),1,0)
  Z3=J(rows(X),1,0)
  Z4=J(rows(X),1,0)
  
  for(i=1; i<=rows(X); i++) {
    p=X[i,1]
	phintWW=X[i,2]
	phintBB=X[i,3]
	theta=X[i,4]
	if ((theta<-0.99)||(theta==4)){
	  V=-99
	}
	else {
      mm_root(V=., &myfunc2(), 0, 5, 0.0001, 1000, p,phintWW,phintBB,theta)
	}
	Z[i]=V
	mm_root(V1=., &myfunc2(), 0, 5, 0.0001, 1000, p,phintWW,phintBB,0.5)
	mm_root(V2=., &myfunc2(), 0, 5, 0.0001, 1000, p,phintWW,phintBB,1)
	mm_root(V3=., &myfunc2(), 0, 5, 0.0001, 1000, p,phintWW,phintBB,1.5)
	mm_root(V4=., &myfunc2(), 0, 5, 0.0001, 1000, p,phintWW,phintBB,2.5)
	Z1[i]=V1
	Z2[i]=V2
	Z3[i]=V3
	Z4[i]=V4
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
end

replace value_ra=. if ((backswitcher==1))
replace value_ra=0 if value_ra<0
replace value_ra=. if theta==-1


*risk-averse indicator:
gen risk_averse=theta>0&backswitcher==0
replace risk_averse=. if backswitcher==1

label var risk_averse "Risk-averse indicator"
label define risk_averse_l 0 "Not risk averse" 1 "Risk-averse"
label values risk_averse risk_averse_l


gen plevel=round(1000*p) //create integer var to code prior probability
gen phintWBs=phintWB*loss
gen phintBWs=phintBW*protectioncost
label var phintWBs "False-neg. prob. x Loss"
label var phintBWs "False-neg. prob. x Prot. cost"




replace bp_val=-bp_val

gen info_effect=bp_val+ip_val //expected difference in earnings between informed and blind protection (subject and decision-specific)
sum info_effect
replace info_effect=0 if info_effect<0 //because wtp is never negative
sum info_effect


gen value_mu=max(0, bp_val-ip_val_mu) //expected diff in earnings between IP and BP accounting for actual subjects' beliefs

hist wtp_time //50% of wtp questions are answered in less than 10.5 seconds

*Calculate discrepancies between WTP and theoretical value for different risk aversion levels:
gen wtp_diff=wtp-value
gen wtp_diff0=wtp-value_ra
gen wtp_diff1=wtp-value_ra1
gen wtp_diff2=wtp-value_ra2
gen wtp_diff3=wtp-value_ra3
gen wtp_diff4=wtp-value_ra4


label values accur_bel accur_bel_l //restoring lost value labels

**Baseline WTP difference: risk-aversion and belief accuracy:
eststo clear
eststo: reg ip_val_diff false_pos false_neg, vce(robust)
eststo: reg ip_val_diff i.plevel false_pos false_neg, vce(robust)
eststo: reg ip_val_diff i.risk_averse##c.false_pos i.risk_averse##c.false_neg, vce(robust)
eststo: reg ip_val_diff i.risk_averse##i.plevel i.risk_averse##c.false_pos i.risk_averse##c.false_neg, vce(robust)
eststo: reg ip_val_diff i.accur_bel##c.false_pos i.accur_bel##c.false_neg, vce(robust)
eststo: reg ip_val_diff i.accur_bel##i.plevel i.accur_bel##c.false_pos i.accur_bel##c.false_neg, vce(robust)
esttab using "./Tables/table_costs1.tex", b(%9.3g) t(%9.1f) ar2(%9.2f) indicate(Prior prob dummies = *.plevel) label title(Expected costs discrepancy) mtitles("" "" "" "" "" "") star("*" 0.10 "**" 0.05 "***" 0.01) nobaselevels compress nogaps replace



eststo clear
eststo: reg ip_val_diff false_pos false_neg if abs(ip_val_diff)<3, vce(robust)
eststo: reg ip_val_diff i.plevel false_pos false_neg if abs(ip_val_diff)<3, vce(robust)
eststo: reg ip_val_diff i.risk_averse##c.false_pos i.risk_averse##c.false_neg if abs(ip_val_diff)<3, vce(robust)
eststo: reg ip_val_diff i.risk_averse##i.plevel i.risk_averse##c.false_pos i.risk_averse##c.false_neg if abs(ip_val_diff)<3, vce(robust)
eststo: reg ip_val_diff i.accur_bel##c.false_pos i.accur_bel##c.false_neg if abs(ip_val_diff)<3, vce(robust)
eststo: reg ip_val_diff i.accur_bel##i.plevel i.risk_averse##c.false_pos i.risk_averse##c.false_neg if abs(ip_val_diff)<3, vce(robust)
esttab using "./Tables/table_costs2.tex", b(%9.3g) t(%9.1f) ar2(%9.2f) indicate(Prior prob dummies = *.plevel) label title(Expected costs discrepancy (without 10\% outliers)) mtitles("" "" "" "" "" "") star("*" 0.10 "**" 0.05 "***" 0.01) nobaselevels compress nogaps replace







**Baseline WTP difference: risk-aversion and belief accuracy:
eststo clear

eststo: reg wtp_diff false_pos false_neg, vce(robust)
eststo: reg wtp_diff i.plevel false_pos false_neg, vce(robust)

eststo: reg wtp_diff i.risk_averse##c.false_pos i.risk_averse##c.false_neg, vce(robust)
eststo: reg wtp_diff i.risk_averse##i.plevel i.risk_averse##c.false_pos i.risk_averse##c.false_neg, vce(robust)

eststo: reg wtp_diff i.accur_bel##c.false_pos i.accur_bel##c.false_neg, vce(robust)
eststo: reg wtp_diff i.accur_bel##i.plevel i.accur_bel##c.false_pos i.accur_bel##c.false_neg, vce(robust)
esttab using "./Tables/table_wtpdiff_01.tex", b(%9.3g) t(%9.1f) ar2(%9.2f) label indicate(Prior dummies=*.plevel) title(WTP for Information (Discrepancy)) mtitles("" "" "" "" "" "" "") star("*" 0.10 "**" 0.05 "***" 0.01) nobaselevels compress nogaps replace




*Demographic variables:
eststo clear
eststo: reg wtp_diff false_pos false_neg, vce(robust)
eststo: reg wtp_diff i.sex##c.false_pos i.sex##c.false_neg, vce(robust)
eststo: reg wtp_diff i.sex##i.plevel i.sex##c.false_pos i.sex##c.false_neg, vce(robust)

eststo: reg wtp_diff i.stat_educ##c.false_pos i.stat_educ##c.false_neg, vce(robust)
eststo: reg wtp_diff i.stat_educ##i.plevel i.stat_educ##c.false_pos i.stat_educ##c.false_neg, vce(robust)

eststo: reg wtp_diff i.old##c.false_pos i.old##c.false_neg, vce(robust)
eststo: reg wtp_diff i.old##i.plevel i.old##c.false_pos i.old##c.false_neg, vce(robust)

eststo: reg wtp_diff i.goodquiz##c.false_pos i.goodquiz##c.false_neg, vce(robust)
eststo: reg wtp_diff i.goodquiz##i.plevel i.goodquiz##c.false_pos i.goodquiz##c.false_neg, vce(robust)


esttab using "./Tables/table_wtpdiff_02.tex", b(%9.3g) t(%9.1f) ar2(%9.2f) label indicate(Prior dummies=*.plevel) title(WTP for Information (Discrepancy, demographic variables)) mtitles("" "" "" "" "" "" "") star("*" 0.10 "**" 0.05 "***" 0.01) nobaselevels compress nogaps replace



*Accounting for risk aversion:
eststo clear
eststo: reg wtp_diff i.plevel false_pos false_neg if !missing(wtp_diff0), vce(robust)
eststo: reg wtp_diff1 i.plevel false_pos false_neg if !missing(wtp_diff0), vce(robust)
eststo: reg wtp_diff2 i.plevel false_pos false_neg if !missing(wtp_diff0), vce(robust)
eststo: reg wtp_diff3 i.plevel false_pos false_neg if !missing(wtp_diff0), vce(robust)
eststo: reg wtp_diff4 i.plevel false_pos false_neg if !missing(wtp_diff0), vce(robust)
eststo: reg wtp_diff0 i.plevel false_pos false_neg if !missing(wtp_diff0), vce(robust)
esttab using "./Tables/table_wtp_ra.tex", b(%9.3g) t(%9.1f) ar2(%9.2f) indicate(Prior dummies=*.plevel) label title(WTP for Information (different risk aversion)) mtitles("$\theta=0$" "$\theta=0.5$" "$\theta=1.0$" "$\theta=1.5$" "$\theta=2.5$" "Heterogeneous $\theta$") star("*" 0.10 "**" 0.05 "***" 0.01) compress nogaps replace




eststo clear
eststo: reg wtp_diff i.plevel false_pos false_neg, vce(robust)
eststo: reg wtp_diff1 i.plevel false_pos false_neg, vce(robust)
eststo: reg wtp_diff2 i.plevel false_pos false_neg, vce(robust)
eststo: reg wtp_diff3 i.plevel false_pos false_neg, vce(robust)
eststo: reg wtp_diff4 i.plevel false_pos false_neg, vce(robust)
esttab using "./Tables/table_wtp_ra2.tex", b(%9.3g) t(%9.1f) ar2(%9.2f) indicate(Prior dummies=*.plevel) label title(WTP for Information (different risk aversion)) mtitles("$\theta=0$" "$\theta=0.5$" "$\theta=1.0$" "$\theta=1.5$" "$\theta=2.5$" "Heterogeneous $\theta$") star("*" 0.10 "**" 0.05 "***" 0.01) compress nogaps replace





*Tobit WTP tables (currently dropped for brevity)





log close
cmdlog using "./Temp/mainwaves_analysis.log"

