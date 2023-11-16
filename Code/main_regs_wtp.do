**********************************************
****-- Main Regressions: WTP FOR INFORMATION --****
**********************************************
set more off
clear all

*!!put your working folder below:
cd C:\Tornado_warnings\Experiment\Alerts_Experiment
*cd C:\Tornado_warnings\Optimal_Alerts

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

hist wtp, title("Distribution of WTP") xtitle("USD") fraction color(navy)
graph export "./Graphs/hist_WTP.png", width(1200) height(800) replace


hist value, title("Distribution of signal's value for a risk-neutral subject") xtitle("USD") fraction note("Value=expected change in costs from BP to IP") color(navy)
graph export "./Graphs/hist_value.png", width(1200) height(800) replace


hist wtp_diff, title("Distribution of WTP discrepancies (WTP - Value)") xtitle("USD") fraction note("Difference between stated wtp and theoretical value for a risk-neutral subject (each choice=obs)") color(navy)
graph export "./Graphs/hist_WTP_discr1.png", width(1200) height(800) replace



*Plotting discrepancies by signal type:
hist wtp_diff if fp_env==0&fn_env==0, title("Distribution of WTP discrepancies (WTP - Value)") xtitle("USD") fraction note("Honest signals") color(navy)
graph export "./Graphs/hist_WTP_discr1h.png", width(1200) height(800) replace


hist wtp_diff if fp_env==1&fn_env==0, title("Distribution of WTP discrepancies (WTP - Value)") xtitle("USD") fraction note("False-positive only signals") color(navy)
graph export "./Graphs/hist_WTP_discr1fp.png", width(1200) height(800) replace

hist wtp_diff if fp_env==0&fn_env==1, title("Distribution of WTP discrepancies (WTP - Value)") xtitle("USD") fraction note("False-negative only signals") color(navy)
graph export "./Graphs/hist_WTP_discr1fn.png", width(1200) height(800) replace

hist wtp_diff if fp_env==1&fn_env==1, title("Distribution of WTP discrepancies (WTP - Value)") xtitle("USD") fraction note("Positive false-positive and false-negative rates") color(navy)
graph export "./Graphs/hist_WTP_discr1fpfn.png", width(1200) height(800) replace



sort subject_id
gen wtp_diff_abs=abs(wtp_diff)
by subject_id: egen totwtp_diff=sum(wtp_diff_abs)
replace totwtp_diff=(1/6)*totwtp_diff
hist totwtp_diff, title("Distribution of WTP discrepancies (WTP - Value)") xtitle("USD") fraction note("Average absolute deviation between stated wtp and theoretical value for a risk-neutral subject, by subject") color(navy)
graph export "./Graphs/hist_WTP_discr2.png", width(1200) height(800) replace


gen highprob=p>0.201
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



**Baseline WTP regs: tobit
gen class_alt=2-class

tab class class_alt
label var class_alt "IP strategy class"
label define clsnames 0 "Bayesians" 1 "Simpletons"
label values class_alt clsnames


gen resp_freqBW=ip_b*freqBW
gen resp_freqBB=ip_b*freqBB
gen resp_freqWB=ip_w*freqWB

gen resp_freqBW1=(1-ip_b)*freqBW
gen resp_freqBB1=(1-ip_b)*freqBB
gen resp_freqWB1=(1-ip_w)*freqWB

*Analyzing the effects of risk aversion on sensitivity of wtp to priors
eststo clear
eststo: reg wtp_diff i.highprob##c.false_neg i.highprob##c.false_pos, vce(cluster subject_id)
eststo: reg wtp_diff i.highprob##c.false_neg i.risk_pref#i.highprob#c.false_neg i.highprob##c.false_pos i.risk_pref#i.highprob#c.false_pos, vce(cluster subject_id)
eststo: reg wtp_diff i.risk_pref##i.highprob##c.false_neg i.risk_pref##i.highprob##c.false_pos, vce(cluster subject_id)
*eststo: reg wtp_diff i.highprob##c.false_neg c.totprot#i.highprob#c.false_neg i.highprob##c.false_pos c.totprot#i.highprob#c.false_pos, vce(cluster subject_id)

eststo: reghdfe wtp_diff i.highprob##c.false_neg i.risk_pref#i.highprob#c.false_neg i.highprob##c.false_pos i.risk_pref#i.highprob#c.false_pos, abs(subject_id) vce(cluster subject_id)
eststo: reghdfe wtp_diff i.risk_pref##i.highprob##c.false_neg i.risk_pref##i.highprob##c.false_pos, abs(subject_id) vce(cluster subject_id)
*eststo: reghdfe wtp_diff i.highprob##c.false_neg c.totprot#i.highprob#c.false_neg i.highprob##c.false_pos c.totprot#i.highprob#c.false_pos, abs(subject_id) vce(cluster subject_id)

esttab using "./Tables/wtp_het_risk.tex", b(%9.3g) se(%9.1f) ar2(%9.2f) label drop(_cons *.risk_pref#*.highprob *.risk_pref#c.false_pos *.risk_pref#c.false_neg) indicate("Full risk pref interactions=*.risk_pref") title(WTP minus Value of Information, risk aversion and sensitivity to FP and FN costs) mtitles("" "" "" "FE" "FE") star("*" 0.10 "**" 0.05 "***" 0.01) nobaselevels compress nogaps replace
esttab using "./Tables/wtp_het_risk_pres.tex", b(%9.3g) ar2(%9.2f) not label drop(_cons *.risk_pref#*.highprob *.risk_pref#*.highprob  *.risk_pref#c.false_pos *.risk_pref#c.false_neg) indicate("Full risk pref interactions=*.risk_pref") mtitles("" "" "" "FE" "FE") star("*" 0.10 "**" 0.05 "***" 0.01) nobaselevels compress nogaps replace

sort subject_id plevel

*Analyzing changes in WTP across prior by subject to find if heterogeneous risk preferences can explain the pattern
save  "./Temp/wtp_discrepancy0.dta", replace

use "./Temp/wtp_discrepancy0.dta", replace
tab plevel phintBW
sort subject_id plevel round
keep subject_id plevel round wtp phintWB phintBW plevel value
reshape wide wtp phintWB phintBW plevel value, i(subject_id) j(round)

gen delta_b11=wtp1-wtp2
gen delta_b12=wtp4-wtp5

gen delta_b21=wtp1-wtp3
gen delta_b22=wtp4-wtp6

gen delta_v11=value1-value2
gen delta_v12=value4-value5

gen delta_v21=value1-value3
gen delta_v22=value4-value6

sum delta_b11 delta_b12
sum delta_b21 delta_b22

ttest delta_b11==delta_b12
ttest delta_b21==delta_b22

gen fp_env2=phintBW2>0
gen fn_env2=phintWB2>0

gen fp_env3=phintBW3>0
gen fn_env3=phintWB3>0
tab fp_env2 fn_env2

ttest delta_b11==delta_b12 if fn_env2==1
ttest delta_v11==delta_v12 if fn_env2==1

ttest delta_b11==delta_b12 if fp_env2==1
ttest delta_v11==delta_v12 if fp_env2==1

ttest delta_b21==delta_b22 if fn_env3==1
ttest delta_v21==delta_v22 if fn_env3==1

ttest delta_b11==delta_b12 if fp_env2==1
ttest delta_v11==delta_v12 if fp_env2==1



use "./Temp/wtp_discrepancy0.dta", replace

gen fp_cost=phintBW*protectioncost
gen fn_cost=phintWB*loss

gen phintWBcat=round(100*phintWB)
gen phintBWcat=round(100*phintBW)

gen false_cost=fp_cost+fn_cost

reg wtp i.plevel false_cost
reg wtp i.plevel fp_cost fn_cost
reg wtp i.plevel i.risk_pref##c.fp_cost i.risk_pref##c.fn_cost


gen delta_p=p-0.275
gen false_pos2=delta_p*fp_cost
gen false_neg2=delta_p*fn_cost

label var fp_cost "FP rate"
label var fn_cost "FN rate"
label var false_pos2 "Prior_change$\times\$FP rate"
label var false_neg2 "Prior_change$\times\$FN rate"


*Testing if subjects' sensitivity to FP/FN costs changes with priors:
eststo clear
eststo: reg wtp i.plevel fp_cost fn_cost false_pos2 false_neg2, vce(cluster subject_id)
estadd scalar adj_r= e(r2_a)
test false_pos2 false_neg2
local p=r(p)
estadd scalar p = `p'
eststo: reg wtp i.plevel fp_cost fn_cost i.plevel#c.fp_cost i.plevel#c.fn_cost, vce(cluster subject_id)
estadd scalar adj_r= e(r2_a)
test 200.plevel#c.fp_cost 200.plevel#c.fn_cost 300.plevel#c.fp_cost 300.plevel#c.fn_cost 500.plevel#c.fp_cost 500.plevel#c.fn_cost
local p=r(p)
estadd scalar p = `p'
eststo: reghdfe wtp i.plevel fp_cost fn_cost false_pos2 false_neg2, absorb(subject_id) vce(cluster subject_id)
estadd scalar adj_r= e(r2_a)
test false_pos2 false_neg2
local p=r(p)
estadd scalar p = `p'
eststo: reghdfe wtp i.plevel fp_cost fn_cost i.plevel#c.fp_cost i.plevel#c.fn_cost, absorb(subject_id) vce(cluster subject_id)
estadd scalar adj_r= e(r2_a)
test 200.plevel#c.fp_cost 200.plevel#c.fn_cost 300.plevel#c.fp_cost 300.plevel#c.fn_cost 500.plevel#c.fp_cost 500.plevel#c.fn_cost
local p=r(p)
estadd scalar p = `p'
esttab using "./Tables/wtp_het_test.tex", b(%9.3g) se(%9.1f) ar2(%9.2f) label drop(_cons fp_cost fn_cost *.plevel) stats(p adj_r N, labels("P(beta_X=0)" "Adjusted \(R^{2}\)" "Observations")) title(WTP for Information (OLS): analyze FP/FN sensitivity by prior) mtitles("" "" "" "FE" "FE" "FE") star("*" 0.10 "**" 0.05 "***" 0.01) nobaselevels compress nogaps replace
esttab using "./Tables/wtp_het_test_pres.tex", b(%9.3g) se(%9.1f) ar2(%9.2f) label drop(_cons fp_cost fn_cost *.plevel) stats(p adj_r N, labels("P(beta_X=0)" "Adjusted \(R^{2}\)" "Observations")) addnotes("Controlling for FP/FN rates, prior, constant omitted") mtitles("" "" "" "FE" "FE" "FE") star("*" 0.10 "**" 0.05 "***" 0.01) nobaselevels compress nogaps replace


*Testing if there is distinguishing btw FP and FN rates:
eststo clear
eststo: reg wtp i.plevel phintWB phintBW, vce(cluster subject_id)
test phintWB=phintBW
estadd scalar adj_r= e(r2_a)
estadd scalar p = r(p)
eststo: reg wtp i.plevel i.plevel#c.phintWB i.plevel#c.phintBW, vce(cluster subject_id)
test  (100.plevel#c.phintWB=100.plevel#c.phintBW) (200.plevel#c.phintWB=200.plevel#c.phintBW)  (300.plevel#c.phintWB=300.plevel#c.phintBW)   (500.plevel#c.phintWB=500.plevel#c.phintBW)
estadd scalar adj_r= e(r2_a)
estadd scalar p = r(p)
eststo: reghdfe wtp i.plevel phintWB phintBW,  absorb(subject_id) vce(cluster subject_id)
test phintWB=phintBW
estadd scalar adj_r= e(r2_a)
estadd scalar p = r(p)
eststo: reghdfe wtp i.plevel i.plevel#c.phintWB i.plevel#c.phintBW,  absorb(subject_id) vce(cluster subject_id)
test  (100.plevel#c.phintWB=100.plevel#c.phintBW) (200.plevel#c.phintWB=200.plevel#c.phintBW)  (300.plevel#c.phintWB=300.plevel#c.phintBW)   (500.plevel#c.phintWB=500.plevel#c.phintBW)
estadd scalar adj_r= e(r2_a)
estadd scalar p = r(p)
esttab using "./Tables/wtp_fpvsfn_test.tex", b(%9.3g) se(%9.1f) ar2(%9.2f) label drop(_cons *.plevel) stats(p adj_r N, labels("P(FP rate=FN rate)" "Adjusted \(R^{2}\)" "Observations")) title(WTP for Information (OLS): comparing sensitivity to FP/FN rates) mtitles("" "" "FE" "FE") star("*" 0.10 "**" 0.05 "***" 0.01) nobaselevels compress nogaps replace
esttab using "./Tables/wtp_fpvsfn_test_pres.tex", b(%9.3g) se(%9.1f) ar2(%9.2f) label drop(_cons *.plevel) stats(p adj_r N, labels("P(FP rate=FN rate)" "Adjusted \(R^{2}\)" "Observations")) mtitles("" "" "FE" "FE") star("*" 0.10 "**" 0.05 "***" 0.01) nobaselevels compress nogaps replace


gen false_prob=phintWB+phintBW
gen false_tot=false_pos+false_neg


**Baseline WTP difference: risk-aversion and belief accuracy:
eststo clear
eststo: reg wtp_diff false_pos false_neg, vce(cluster subject_id)
eststo: reghdfe wtp_diff false_pos false_neg, abs(subject_id) vce(cluster subject_id)
eststo: reghdfe wtp_diff i.risk_pref##c.false_pos i.risk_pref##c.false_neg, abs(subject_id) vce(cluster subject_id)
eststo: reghdfe wtp_diff i.accur_bel2##c.false_pos i.accur_bel2##c.false_neg, abs(subject_id) vce(cluster subject_id)
eststo: reghdfe wtp_diff i.plevel##c.false_pos i.plevel##c.false_neg, abs(subject_id) vce(cluster subject_id)
*eststo: reghdfe wtp_diff i.accur_bel2#i.risk_pref##c.false_pos i.accur_bel2#i.risk_pref##c.false_neg if plevel<300, abs(subject_id) vce(cluster subject_id)
esttab using "./Tables/table_wtpdiff_01rfull_ag.tex", b(%9.3g) se(%9.1f) ar2(%9.2f) label title(WTP minus Value of Information (OLS)) mtitles("" "" "" "" "") star("*" 0.10 "**" 0.05 "***" 0.01) nobaselevels compress nogaps replace
esttab using "./Tables/table_wtpdiff_01rfull_ag_pres.tex", b(%9.3g) ar2(%9.2f) not indicate("Risk pref dummies = *.risk_pref" "Prior dummies = *.plevel" "Belief accuracy = *.accur_bel2" ) label mtitles("" "FE" "FE" "FE" "FE") star("*" 0.10 "**" 0.05 "***" 0.01) nobaselevels compress nogaps replace
esttab using "./Tables/table_wtpdiff_01r_ag.tex", b(%9.3g) se(%9.1f) ar2(%9.2f) label title(WTP minus Value of Information (OLS)) mtitles("" "" "" "" "") star("*" 0.10 "**" 0.05 "***" 0.01) noomit nobaselevels compress nogaps replace





*Demographic variables:
eststo clear
*eststo: reg wtp_diff false_pos false_neg, vce(cluster subject_id)
eststo: reg wtp_diff false_pos false_neg i. sex i.sex#c.false_pos i.sex#c.false_neg, vce(cluster subject_id)
eststo: reg wtp_diff i.sex##i.plevel i.sex##c.false_pos i.sex##c.false_neg, vce(cluster subject_id)
eststo: reg wtp_diff i.stat_educ##c.false_pos i.stat_educ##c.false_neg, vce(cluster subject_id)
eststo: reg wtp_diff i.stat_educ##i.plevel i.stat_educ##c.false_pos i.stat_educ##c.false_neg, vce(cluster subject_id)
eststo: reg wtp_diff i.old##c.false_pos i.old##c.false_neg, vce(cluster subject_id)
eststo: reg wtp_diff i.old##i.plevel i.old##c.false_pos i.old##c.false_neg, vce(cluster subject_id)
esttab using "./Tables/table_wtpdiff_02.tex", b(%9.3g) se(%9.1f) ar2(%9.2f) label indicate(Prior dummies=*.plevel) title(WTP minus Value of Information: demographic determinants) mtitles("" "" "" "" "" "" "" "" "") star("*" 0.10 "**" 0.05 "***" 0.01) nobaselevels compress nogaps replace


*Prior heterogeneity and instructions understanding:
eststo clear
eststo: reg wtp_diff i.plevel i.highprob##c.false_neg i.goodquiz#i.highprob#c.false_neg i.highprob##c.false_pos i.goodquiz#i.highprob#c.false_pos, vce(cluster subject_id)
eststo: reg wtp_diff i.plevel i.goodquiz##i.highprob##c.false_neg i.goodquiz##i.highprob##c.false_pos, vce(cluster subject_id)
*eststo: reg wtp_diff i.highprob##c.false_neg c.totprot#i.highprob#c.false_neg i.highprob##c.false_pos c.totprot#i.highprob#c.false_pos, vce(cluster subject_id)

eststo: reghdfe wtp_diff i.plevel i.highprob##c.false_neg i.goodquiz#i.highprob#c.false_neg i.highprob##c.false_pos i.goodquiz#i.highprob#c.false_pos, abs(subject_id) vce(cluster subject_id)
eststo: reghdfe wtp_diff i.plevel i.goodquiz##i.highprob##c.false_neg i.goodquiz##i.highprob##c.false_pos, abs(subject_id) vce(cluster subject_id)

esttab using "./Tables/table_wtpdiff_instr.tex", b(%9.3g) se(%9.1f) ar2(%9.2f) label indicate(Prior dummies=*.plevel) title(WTP minus Value of Information: instructions understanding) mtitles("" "" "" "" "" "" "" "" "") star("*" 0.10 "**" 0.05 "***" 0.01) nobaselevels compress nogaps replace



eststo clear
eststo: reg wtp_diff i.highprob##c.false_neg i.highprob##c.false_pos if goodquiz==0, vce(cluster subject_id)
eststo: reg wtp_diff i.highprob##c.false_neg i.highprob##c.false_pos if goodquiz==1, vce(cluster subject_id)

eststo: reg wtp_diff i.plevel i.highprob##c.false_neg i.highprob##c.false_pos if goodquiz==0, vce(cluster subject_id)
eststo: reg wtp_diff i.plevel i.highprob##c.false_neg i.highprob##c.false_pos if goodquiz==1, vce(cluster subject_id)
*eststo: reg wtp_diff i.highprob##c.false_neg c.totprot#i.highprob#c.false_neg i.highprob##c.false_pos c.totprot#i.highprob#c.false_pos, vce(cluster subject_id)
esttab using "./Tables/table_wtpdiff_instr.tex", b(%9.3g) se(%9.1f) ar2(%9.2f) label indicate(Prior dummies=*.plevel) title(WTP minus Value of Information: instructions understanding) mtitles("Bad" "Good" "Bad" "Good") star("*" 0.10 "**" 0.05 "***" 0.01) nobaselevels compress nogaps replace



*By prior prob for FP and FN costs:
eststo clear
eststo: reghdfe wtp_diff false_pos false_neg if plevel==100, abs(subject_id) vce(cluster subject_id)
eststo: reghdfe wtp_diff false_pos false_neg if plevel==200, abs(subject_id) vce(cluster subject_id)
eststo: reghdfe wtp_diff false_pos false_neg if plevel==300, abs(subject_id) vce(cluster subject_id)
eststo: reghdfe wtp_diff false_pos false_neg if plevel==500, abs(subject_id) vce(cluster subject_id)
esttab using "./Tables/table_wtpdiff_03.tex", b(%9.3g) se(%9.1f) ar2(%9.2f) label title(WTP - Value of Information, by prior) addnotes("Subject fixed effects are included.") mtitles("0.1" "0.2" "0.3" "0.5") star("*" 0.10 "**" 0.05 "***" 0.01) nobaselevels compress nogaps replace


*By prior prob for FP and FN rates:
eststo clear
eststo: tobit wtp phintBW phintWB if plevel==100, ll(0) ul(5)
test phintBW=phintWB
local p=r(p)
estadd scalar p = `p'
estadd scalar adj_r= e(r2_a)
eststo: tobit wtp phintBW phintWB if plevel==200, ll(0) ul(5)
test phintBW=phintWB
local p=r(p)
estadd scalar p = `p'
estadd scalar adj_r= e(r2_a)
eststo: tobit wtp phintBW phintWB if plevel==300, ll(0) ul(5)
test phintBW=phintWB
local p=r(p)
estadd scalar p = `p'
estadd scalar adj_r= e(r2_a)
eststo: tobit wtp phintBW phintWB if plevel==500, ll(0) ul(5)
test phintBW=phintWB
local p=r(p)
estadd scalar p = `p'
estadd scalar adj_r= e(r2_a)
esttab using "./Tables/table_wtpdiff_04tob.tex", b(%9.3g) se(%9.1f) ar2(%9.2f) label stats(p adj_r N, labels("P(FP rate=FN rate)" "Adjusted \(R^{2}\)" "Observations")) title(WTP for Information, by prior (tobit)) mtitles("0.1" "0.2" "0.3" "0.5") star("*" 0.10 "**" 0.05 "***" 0.01) nobaselevels compress nogaps replace
esttab using "./Tables/table_wtpdiff_04tob_pres.tex", b(%9.3g) se(%9.1f) ar2(%9.2f) drop(_cons) stats(p adj_r N, labels("P(FP rate=FN rate)" "Adjusted \(R^{2}\)" "Observations")) eqlabels(none) label mtitles("0.1" "0.2" "0.3" "0.5") addnotes("Tobit regression of WTP, constant omitted") star("*" 0.10 "**" 0.05 "***" 0.01) nobaselevels compress nogaps replace




save "./Temp/wtp_discrepancy0.dta", replace

use "./Temp/wtp_discrepancy0.dta", replace
drop highprob
gen highprob=p>0.201


**Calculate average WTP discrepancy by signal type
tempname p1
postfile `p1' false_pos false_neg wtp_diff ptest using "./Temp/wtp_by_environment.dta", replace


forvalues i=0/1{
  forvalues j=0/1{
	    ttest wtp_diff==0 if fp_env==`i'&fn_env==`j', level(95)
		local ptest=r(p)
		local wtp_diff=r(mu_1)
		post `p1' (`i') (`j') (`wtp_diff') (`ptest')
    }
  }

postclose `p1'


**Calculate average WTP discrepancy by signal type
tempname p2
postfile `p2' highprob false_pos false_neg wtp_diff ptest using "./Temp/wtp_by_environment_det.dta", replace

forvalues k=0/1{
  forvalues i=0/1{
    forvalues j=0/1{
	    ttest wtp_diff==0 if highprob==`k'&fp_env==`i'&fn_env==`j', level(95)
		local ptest=r(p)
		local wtp_diff=r(mu_1)
		post `p2' (`k') (`i') (`j') (`wtp_diff') (`ptest')
    }
  }
}
postclose `p2'

use "./Temp/wtp_by_environment.dta", replace

tostring false_pos false_neg, replace
foreach var of varlist false_pos false_neg{
  replace `var'="No" if `var'=="0"
  replace `var'="Yes" if `var'=="1"
}
bro
format wtp_diff ptest %9.3f
listtex using "./Tables/bigpicture_wtp.tex", type rstyle(tabular) head("\begin{table}[H]\centering \caption{Average WTP discrepancy (WTP-Value) by Signal Type} \begin{tabular}{cccc} \hline \hline" `"\textbf{False-positive}&\textbf{False-negative}&\textbf{Mean WTP discrepancy}& \textbf{P($=0$)}\\ \hline"') foot("\hline \end{tabular} \end{table}") replace
listtex using "./Tables/bigpicture_wtp_pres.tex", type rstyle(tabular) head("\begin{table}[H]\centering \begin{tabular}{cccc} \hline \hline" `"\textbf{False-positive}&\textbf{False-negative}&\textbf{Mean WTP discrepancy}& \textbf{P($=0$)}\\ \hline"') foot("\hline \end{tabular} \end{table}") replace

use "./Temp/wtp_by_environment_det.dta", replace

tostring false_pos false_neg, replace
foreach var of varlist false_pos false_neg highprob{
  replace `var'="No" if `var'=="0"
  replace `var'="Yes" if `var'=="1"
}
bro
format wtp_diff ptest %9.3f
listtex using "./Tables/bigpicture_wtp_det.tex", type rstyle(tabular) head("\begin{table}[H]\centering \caption{Average WTP discrepancy (WTP-Value) by Signal Type} \begin{tabular}{cccc} \hline \hline" `"\textbf{False-positive}&\textbf{False-negative}&\textbf{Mean WTP discrepancy}& \textbf{P($=0$)}\\ \hline"') foot("\hline \end{tabular} \end{table}") replace
listtex using "./Tables/bigpicture_wtp_det_pres.tex", type rstyle(tabular) head("\begin{table}[H]\centering \begin{tabular}{cccc} \hline \hline" `"\textbf{False-positive}&\textbf{False-negative}&\textbf{Mean WTP discrepancy}& \textbf{P($=0$)}\\ \hline"') foot("\hline \end{tabular} \end{table}") replace




use "./Temp/wtp_discrepancy0.dta", replace
*How many switch choices after changing prior probabilities:
gen wtp_p_change=wtp~=wtp[_n-3]
replace wtp_p_change=0 if round<4
bys subject_id: egen nchanges=sum(wtp_p_change)

gen wtp_s_change=wtp~=wtp[_n-1]
replace wtp_s_change=0 if inlist(round,1,2,4)
bys subject_id: egen nchanges_s=sum(wtp_s_change)
tab nchanges_s

**Removing subjects which do not update their WTP after increasing priors:
eststo clear
eststo: reghdfe wtp_diff phintBW phintWB if nchanges>0, abs(subject_id) vce(cluster subject_id)
eststo: reghdfe wtp_diff phintBW phintWB if plevel==100&nchanges>0, abs(subject_id) vce(cluster subject_id)
eststo: reghdfe wtp_diff phintBW phintWB if plevel==200&nchanges>0, abs(subject_id) vce(cluster subject_id)
eststo: reghdfe wtp_diff phintBW phintWB if plevel==300&nchanges>0, abs(subject_id) vce(cluster subject_id)
eststo: reghdfe wtp_diff phintBW phintWB if plevel==500&nchanges>0, abs(subject_id) vce(cluster subject_id)
esttab using "./Tables/table_wtpdiff_06s.tex", b(%9.3g) se(%9.1f) ar2(%9.2f) label title(WTP - Value of Information, by prior) mtitles("All" "0.1" "0.2" "0.3" "0.5")  addnotes("Only subjects who change their decisions across priors") star("*" 0.10 "**" 0.05 "***" 0.01) nobaselevels compress nogaps replace



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

