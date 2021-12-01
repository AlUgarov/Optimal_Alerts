
**ANALYZE THE FIRST WAVE (planned 100 participants)***
set more off
clear all

**Requires: estout, moremata

*!!put your working folder below:
cd C:\Tornado_warnings\Optimal_Alerts


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
import delimited using "./Input/Data_All_mainwaves.csv", rowrange(4:109) varnames(nonames) clear
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

*identify correct quiz answers
gen blind_correct=(q157==2)+(q158==3)
tab blind_correct
gen informed_correct=(q120==2)+(q119==2)+(q121==1)+(q118==3)+(q117==4)
tab informed_correct
gen wtpq_correct=(q164==3)+(q166==3)
tab wtpq_correct
gen ncorrect=blind_correct+informed_correct+wtpq_correct
tab ncorrect

encode participant_id, gen(subject_id) //as participant_id is initially a string

save "./Temp/mainwaves_wide.dta", replace


***SANITY CHECKS*******************

**??**



**BLIND PROTECTION ANALYSIS**
* bp - protection decision (0 - do not protect, 1 - protect)
keep subject_id bp_* bp_time_*
reshape long bp_ bp_time_,  i(subject_id) j(round)
rename bp_ bp
rename bp_time_ submittime
gen p=0.05*round
reg bp p

sum submittime

**Identify switchers: change to no protection for higher probabilities
sort subject_id round
gen backswitch=(bp==0)&(bp[_n-1]==1)
replace backswitch=0 if round==1
by subject_id: egen backswitcher=max(backswitch)
by subject_id: egen nbswitches=sum(backswitch)
save "./Temp/blind_prot_temp.dta", replace


gen backswitchround=round*backswitch

replace backswitchround=. if backswitchround==0

gen loss=20
gen protectioncost=5
gen bp_val=-((1-bp)*p*loss+bp*protectioncost) //expected costs of blind protection decision
sum bp_val
sort subject_id round
list subject_id round bp if backswitcher==1


*repairing back switchers
*if only one back switch and it can be repaired by a single change
*then change the switchround
*repair 1: switching to no protection in the last round
*repair 2: switching to protection 
gen bpc=bp
replace bpc=1-bpc if (nbswitches==1)&(bpc[_n-1]!=bpc)&(bpc[_n+1]!=bpc)&(round>1)&(round<6)
replace bpc=1-bpc if (nbswitches==1)&(bpc[_n-1]!=bpc)&(bpc==0)&(round==6)
list subject_id round bp bpc if backswitcher==1

gen switch=(bpc==1)&(bpc[_n-1]==0)
replace switch=0 if round==1
by subject_id: egen switcher=max(switch)
gen switchround=round*switch
replace switchround=. if switchround==0

drop backswitch backswitcher
gen backswitch=(bpc==0)&(bpc[_n-1]==1)
replace backswitch=0 if round==1
by subject_id: egen backswitcher=max(backswitch)


save "./Temp/bp_val.dta", replace


*Collapsing to have one obs per participant (study participant's characteristics):
collapse (mean) bp submittime (sum) totprot=bp (max) switcher backswitcher nbswitches (min) maxspeed=submittime firstswitch=switchround backswitchround, by(subject_id)

tab switcher
tab backswitcher
tab firstswitch //the distribution of the first switching round
tab backswitchround

scatter submittime backswitcher


**Estimate subjects' risk aversion from the blind protection choices:
gen switchprob=0.05*firstswitch
replace switchprob=-0.5 if missing(switchprob)
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
  (X, Z)
  mata drop funk()
  idx = st_addvar("float", "theta")
  st_store(., "theta", Z)
end

replace theta=. if backswitcher==1
replace theta=9 if (totprot==6)

//Make the table on distribution of thetas
tab theta


* the oneway table(doesn't work so far)
*eststo clear
*eststo: estpost tabstat theta, by(theta)
*esttab using "./Tables/thetas.tex", cell(b) unstack b(%9.3g) t(%9.1f) ar2(%9.2f) title ("Relative risk avers distribution") label replace

gen byte allprotect=(totprot==6)
save "./Temp/blind_collapsed.dta", replace




**CHANGE TO THE PANEL STRUCTURE FOR THE OTHER TASKS********************
**Panel: participant_id round 
**because all the tasks except the blind protection have the same N of rounds (6) 
use "./Temp/mainwaves_wide.dta", replace

reshape long ip_w_ ip_b_ ip_time_ be_w_ be_b_ be_time_w_ be_time_b_ wtp_1_ wtp_2_ wtp_3_ wtp_4_ wtp_5_ wtp_6_ wtp_7_ wtp_8_ wtp_9_ wtp_10_ wtp_11_ wtp_time_, i(subject_id participant_id) j(round)
rename *_ *
rename ip_time time_ip

***Merging the treatment sequences***
merge m:1 seq round using "./Temp/exp_treatments.dta"
tab _merge
drop _merge

**Merging risk aversion vars from the blind protection task
merge m:1 subject_id using "./Temp/blind_collapsed.dta"
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
replace be_w=0.01*(100-be_w) //rescale belief elicitation responses to probabilities
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

gen ip_val=-(p*(phintWB*(1-ip_w)+phintBB*(1-ip_b))*loss+p*(phintWB*ip_w+phintBB*ip_b)*protectioncost+(1-p)*(phintWW*ip_w+phintBW*ip_b)*protectioncost)
sum ip_val

//Calculate the optimal protection strategy based on cost-loss ratio:
gen ip_w_o=0
replace ip_w_o=1 if post_probW>=protectioncost/loss
gen ip_b_o=0
replace ip_b_o=1 if post_probB>=protectioncost/loss

//Calculate exp costs under the optimal strategy:
gen ip_val_o=-(p*(phintWB*(1-ip_w_o)+phintBB*(1-ip_b_o))*loss+p*(phintWB*ip_w_o+phintBB*ip_b_o)*protectioncost+(1-p)*(phintWW*ip_w_o+phintBW*ip_b_o)*protectioncost)

reg ip_val ip_val_o

label var ip_val "Exp. costs"
label var ip_val_o "Optimal exp. costs"


*Saving the cleaned dataset with the panel structure
save "./Output/main_waves.dta", replace


*Remove the timing information for now
drop *click* *page* history q104-q106 q114 q115
rename post_probW post_probw
rename post_probB post_probb
keep subject_id round ip_w ip_b be_w be_b be_time_w be_time_b post_probw post_probb time_ip honest ncorrect p phintBW phintWB phintBB

*Reshape to (subject_id round hint) long format
reshape long ip_ be_ be_time_ post_prob, i(subject_id round time_ip honest ncorrect p phintBW phintWB phintBB) j(hint) string
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
esttab using "./Tables/table_ip1.tex", b(%9.3g) t(%9.1f) aic(%9.2f) label title(Informed Protection) mtitles("All" "All" "Smart" "Smart") star("*" 0.10 "**" 0.05 "***" 0.01) compress nogaps replace

**response to elicited (and posterior) probabilities of black***
eststo clear
eststo: probit ip_ be_, vce(robust)
eststo: probit ip_ be_ post_prob, vce(robust)
eststo: probit ip_ be_ post_prob if ncorrect>6, vce(robust)
esttab using "./Tables/table_ip2.tex", b(%9.3g) t(%9.1f) aic(%9.2f) label title(Informed Protection: Response to Reported Beliefs) mtitles("All" "All" "Smart") star("*" 0.10 "**" 0.05 "***" 0.01) compress nogaps replace

**testing for timing and correlation
xtreg ip_ be_, fe vce(robust)
xtreg ip_ be_ if time_ip<20, fe vce(robust)
reg ip_ post_prob, vce(robust)
reg ip_ post_prob if time_ip<20, vce(robust)



****-- BELIEF ELICITATION --****

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


*Testing the accuracy of reported beliefs***
eststo clear
eststo: reg be_ post_prob, vce(robust)
eststo: reg be_ post_prob if ncorrect>6, vce(robust)
eststo: reg be_ post_prob if (honest==0), vce(robust)
esttab using "./Tables/table_be1.tex", b(%9.3g) t(%9.1f) ar2(%9.2f) label title(Belief Elicitation: Belief vs Posterior) mtitles("All" "Not_honest" "Good quiz") star("*" 0.10 "**" 0.05 "***" 0.01) compress nogaps replace



*Decomposition of belief updating (coeffs should be all ones)***
eststo clear
eststo: reg lt_bel lt_prior signalB signalW, vce(robust)
eststo: xtreg lt_bel lt_prior signalB signalW, fe vce(robust)
eststo: xtreg lt_bel lt_prior signalB signalW if ncorrect>6, fe vce(robust)
esttab using "./Tables/table_be3.tex", b(%9.3g) t(%9.1f) ar2(%9.2f) label title(Belief Elicitation: Decomposition) mtitles("OLS" "FE" "Smart, FE") star("*" 0.10 "**" 0.05 "***" 0.01) compress nogaps replace



****-- WTP FOR INFORMATION --****
use "./Output/main_waves.dta", replace
xtset subject_id round

*calculate theoretical value of signal for a risk-neutral subject:
gen cost_bp=min(p*loss, protectioncost) //expected blind protection cost
gen false_pos=(p*phintBB+(1-p)*phintBW)*protectioncost //expected protection costs (protect as long as there is a signal)
gen false_neg=p*phintWB*loss //expected false positive costs
gen cost_ip=false_neg+false_pos //total informed protection costs
gen value=max(0, cost_bp-cost_ip) //theoretical value for a risk-neutral subject

label var false_pos "Prot. costs"
label var false_neg "False neg. costs"
label var cost_bp "BP costs"
label var phintWB "False neg. rate"
label var phintBW "False pos. rate"

*uniform CRRA values to calculate the signal's values later
gen theta1=0.5
gen theta2=1
gen theta3=1.5
gen theta4=2.5

*calculate theoretical value for a risk-averse subject (both heterogeneous and uniform CRRA coeffs):
mata:
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

*belief-adjusted value for a risk-neutral subject:
*gen value_b=max(0,-(p*phintWB*loss-(p*phintBB+(1-p)*phintBW)*protectioncost))



/*Merging blind protection choices to get expected costs under blind protection for each probability*/
merge m:1 subject_id p using "./Temp/bp_val.dta"

gen info_effect=ip_val-bp_val //expected difference in earnings between informed and blind protection (subject and decision-specific)
sum info_effect
replace info_effect=0 if info_effect<0 //because wtp is never negative


sum value wtp wtp_time
hist wtp_time //Most people answer wtp question in less than 10 seconds

*Calculate discrepancies between WTP and theoretical value for different risk aversion levels:
gen wtp_diff=wtp-value
gen wtp_diff0=wtp-value_ra
gen wtp_diff1=wtp-value_ra1
gen wtp_diff2=wtp-value_ra2
gen wtp_diff3=wtp-value_ra3
gen wtp_diff4=wtp-value_ra4

*Estimate the discrepancy between wtp and risk-neutral value as a function of false positive losses, expected protection costs and blind protection costs:
eststo clear
eststo: reg wtp_diff false_pos false_neg, vce(robust)
eststo: reg wtp_diff false_pos false_neg if theta>0&backswitcher==0, vce(robust)
eststo: reg wtp_diff false_pos false_neg if theta<0&backswitcher==0, vce(robust)
eststo: reg wtp_diff false_pos false_neg if backswitcher==1, vce(robust)
esttab using "./Tables/table_wtpdiff.tex", b(%9.3g) t(%9.1f) ar2(%9.2f) label title(WTP for Information (Discrepancy)) mtitles("All" "Risk-averse" "Risk-loving" "Switchers") star("*" 0.10 "**" 0.05 "***" 0.01) compress nogaps replace


*Estimate wtp as a function of false positive losses, expected protection costs and blind protection costs:
eststo clear
eststo: tobit wtp cost_bp false_pos false_neg, ll(0) ul(5)
eststo: tobit wtp cost_bp false_pos false_neg if theta>0&backswitcher==0, ll(0) ul(5)
eststo: tobit wtp cost_bp false_pos false_neg if theta<0&backswitcher==0, ll(0) ul(5)
eststo: tobit wtp cost_bp false_pos false_neg if backswitcher==1, ll(0) ul(5)
esttab using "./Tables/table_wtp01.tex", b(%9.3g) t(%9.1f) aic(%9.2f) label title(WTP for Information (Tobit Estimation)) mtitles("All" "Risk-averse" "Risk-loving" "Switchers") star("*" 0.10 "**" 0.05 "***" 0.01) compress nogaps replace




*Compare coefficient estimates between risk-averse and risk-loving subjects
tobit wtp cost_bp false_pos false_neg if theta>0&backswitcher==0, ll(0) ul(5)
est store riskav
tobit wtp cost_bp false_pos false_neg if theta<0&backswitcher==0, ll(0) ul(5)
est store risklov
suest riskav risklov
test [riskav_model]cost_bp=[risklov_model]cost_bp
test [riskav_model]false_pos=[risklov_model]false_pos
test [riskav_model]false_neg=[risklov_model]false_neg
test [riskav_model=risklov_model]



*Accounting for risk aversion:
eststo clear
eststo: reg wtp_diff0 cost_bp false_pos false_neg, vce(robust)
eststo: reg wtp_diff1 cost_bp false_pos false_neg, vce(robust)
eststo: reg wtp_diff2 cost_bp false_pos false_neg, vce(robust)
eststo: reg wtp_diff3 cost_bp false_pos false_neg, vce(robust)
eststo: reg wtp_diff4 cost_bp false_pos false_neg, vce(robust)
esttab using "./Tables/table_wtp_ra.tex", b(%9.3g) t(%9.1f) ar2(%9.2f) label title(WTP for Information (different risk aversion)) mtitles("Heterogeneous" "$\theta=0.5$" "$\theta=1.0$" "$\theta=1.5$" "$\theta=2.5$") star("*" 0.10 "**" 0.05 "***" 0.01) compress nogaps replace

*Actual expected costs vs optimal theoretical costs of informed protection:
eststo clear
eststo: reg ip_val ip_val_o, vce(robust)
eststo: reg  ip_val ip_val_o p, vce(robust)
eststo: reg  ip_val ip_val_o p phintWB phintBW, vce(robust)
eststo: xtreg ip_val ip_val_o, fe vce(robust)
eststo: xtreg  ip_val ip_val_o p, fe vce(robust)
eststo: xtreg  ip_val ip_val_o phintWB phintBW, fe vce(robust)
esttab using "./Tables/table_costs0.tex", b(%9.3g) t(%9.1f) ar2(%9.2f) label title(Actual Exp. Costs vs Theoretical Costs) mtitles("OLS" "OLS" "OLS" "FE" "FE" "FE") star("*" 0.10 "**" 0.05 "***" 0.01) compress nogaps replace


log close
cmdlog using "./Temp/mainwaves_analysis.log"

