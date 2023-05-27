**********************************************
****-- Main Regressions: WTP FOR INFORMATION --****
**********************************************
global ROOT "C:\Users\gaduh\Documents\GitHub\Optimal_Alerts"


clear all
use "$ROOT/Output/main_waves.dta", clear
xtset subject_id round

merge m:1 subject_id round using "$ROOT/Temp/bel_accuracy.dta" //average belief accuracy by subject
	drop _merge


merge m:1 participant_id p using "$ROOT/Temp/bp_val.dta" //risk aversion, blind prot choices and demographic vars
	drop if _merge==2
	drop _merge

merge m:1 subject_id using "$ROOT/Temp/ipclasses.dta" //risk aversion, blind prot choices and demographic vars
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
label value accur_bel3 accur_bel_l


sort subject_id
gen wtp_diff_abs=abs(wtp_diff)
by subject_id: egen totwtp_diff=sum(wtp_diff_abs)
replace totwtp_diff=(1/6)*totwtp_diff

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
label define risk_prefl 0 "Risk-neutral" 1 "Risk-loving" 2 "Risk-averse" 3 "Unmeasured risk"
label value risk_pref risk_prefl

gen accur_bel1=1-accur_bel
label var accur_bel1 "Beliefs accuracy"
label def accur_bel1l 0 "Accurate belief" 1 "Inaccurate belief"
label value accur_bel1 accur_bel1l

gen inac_bel3=1-round(accur_bel3)
label def inac_bel3 0 "Accurate Subject-by-Treatment" 1 "Inaccurate Subject-by-Treatment"
label value inac_bel3 inac_bel3

bys subject_id: egen inac_count=sum(inac_bel3)
sum inac_count, detail
gen inac_sub = inac_count > r(p50)
label def inac_sub3 0 "Accurate Subject" 1 "Inaccurate Subject"
label value inac_sub inac_sub3

*Self-reported strategies
*gen strategy_short=strategy_ip
*replace strategy_short=3 if strategy_short>3
label var strategy_short "Strategy"
label define strategy_short_l 1 "Rational" 2 "Seek honest" 3 "Other"
label values strategy_short strategy_short_l


label var be_change "Belief change"
label var confid "Certainty"

bys fp_env fn_env: ttest wtp_diff == 0 if plevel<300, level(95)

**Baseline WTP difference: risk-aversion and belief accuracy:
eststo clear

/*eststo: reg wtp_diff false_pos false_neg, vce(cluster subject_id)
	test false_pos = false_neg
	estadd scalar test_eq=r(p)
*/
eststo: reghdfe wtp_diff false_pos false_neg, abs(subject_id) vce(cluster subject_id)
	test false_pos = false_neg
		estadd scalar test_eq=r(p)
	
eststo: reghdfe wtp_diff i.risk_pref##c.false_pos i.risk_pref##c.false_neg, abs(subject_id) vce(cluster subject_id)
	test false_pos = false_neg
		estadd scalar test_eq=r(p)
	lincom false_pos+1.risk_pref#c.false_pos
		estadd scalar b_RL_fp_X=r(estimate)
		estadd scalar se_RL_fp_X=r(se)
		estadd scalar p_RL_fp_X=r(p)
	lincom false_neg+1.risk_pref#c.false_neg
		estadd scalar b_RL_fn_X=r(estimate)
		estadd scalar se_RL_fn_X=r(se)
		estadd scalar p_RL_fn_X=r(p)
	lincom false_pos+2.risk_pref#c.false_pos
		estadd scalar b_RA_fp_X=r(estimate)
		estadd scalar se_RA_fp_X=r(se)
		estadd scalar p_RA_fp_X=r(p)
	lincom false_neg+2.risk_pref#c.false_neg
		estadd scalar b_RA_fn_X=r(estimate)
		estadd scalar se_RA_fn_X=r(se)
		estadd scalar p_RA_fn_X=r(p)
	
eststo: reghdfe wtp_diff i.inac_bel3##c.false_pos i.inac_bel3##c.false_neg, abs(subject_id) vce(cluster subject_id)
	test false_pos = false_neg
		estadd scalar test_eq=r(p)

eststo: reghdfe wtp_diff i.inac_sub##c.false_pos i.inac_sub##c.false_neg, abs(subject_id) vce(cluster subject_id)
	test false_pos = false_neg
		estadd scalar test_eq=r(p)

		
eststo: reghdfe wtp_diff i.risk_pref##c.false_pos i.risk_pref##c.false_neg if !inac_sub, abs(subject_id) vce(cluster subject_id)
	test false_pos = false_neg
		estadd scalar test_eq=r(p)
	lincom false_pos+1.risk_pref#c.false_pos
		estadd scalar b_RL_fp_X=r(estimate)
		estadd scalar se_RL_fp_X=r(se)
		estadd scalar p_RL_fp_X=r(p)
	lincom false_neg+1.risk_pref#c.false_neg
		estadd scalar b_RL_fn_X=r(estimate)
		estadd scalar se_RL_fn_X=r(se)
		estadd scalar p_RL_fn_X=r(p)
	lincom false_pos+2.risk_pref#c.false_pos
		estadd scalar b_RA_fp_X=r(estimate)
		estadd scalar se_RA_fp_X=r(se)
		estadd scalar p_RA_fp_X=r(p)
	lincom false_neg+2.risk_pref#c.false_neg
		estadd scalar b_RA_fn_X=r(estimate)
		estadd scalar se_RA_fn_X=r(se)
		estadd scalar p_RA_fn_X=r(p)
		
eststo: reghdfe wtp_diff i.risk_pref##c.false_pos i.risk_pref##c.false_neg if inac_sub, abs(subject_id) vce(cluster subject_id)
	test false_pos = false_neg
		estadd scalar test_eq=r(p)
	lincom false_pos+1.risk_pref#c.false_pos
		estadd scalar b_RL_fp_X=r(estimate)
		estadd scalar se_RL_fp_X=r(se)
		estadd scalar p_RL_fp_X=r(p)
	lincom false_neg+1.risk_pref#c.false_neg
		estadd scalar b_RL_fn_X=r(estimate)
		estadd scalar se_RL_fn_X=r(se)
		estadd scalar p_RL_fn_X=r(p)
	lincom false_pos+2.risk_pref#c.false_pos
		estadd scalar b_RA_fp_X=r(estimate)
		estadd scalar se_RA_fp_X=r(se)
		estadd scalar p_RA_fp_X=r(p)
	lincom false_neg+2.risk_pref#c.false_neg
		estadd scalar b_RA_fn_X=r(estimate)
		estadd scalar se_RA_fn_X=r(se)
		estadd scalar p_RA_fn_X=r(p)

		
/*		
eststo: reghdfe wtp_diff i.inac_bel3##i.risk_pref##c.false_pos i.inac_bel3##i.risk_pref##c.false_neg, abs(subject_id) vce(cluster subject_id)
	test false_pos = false_neg
		estadd scalar test_eq=r(p)
	lincom false_pos+1.risk_pref#c.false_pos
		estadd scalar b_RL_fp_AB=r(estimate)
		estadd scalar se_RL_fp_AB=r(se)
	lincom false_neg+1.risk_pref#c.false_neg
		estadd scalar b_RL_fn_AB=r(estimate)
		estadd scalar se_RL_fn_AB=r(se)
	lincom false_pos+2.risk_pref#c.false_pos
		estadd scalar b_RA_fp_AB=r(estimate)
		estadd scalar se_RA_fp_AB=r(se)
	lincom false_neg+2.risk_pref#c.false_neg
		estadd scalar b_RA_fn_AB=r(estimate)
		estadd scalar se_RA_fn_AB=r(se)

*/		
/*eststo: reghdfe wtp_diff i.plevel##c.false_pos i.plevel##c.false_neg, abs(subject_id) vce(cluster subject_id)
	test false_pos = false_neg
	estadd scalar test_eq=r(p)
*/

/*
#delimit ;
	esttab * using "$ROOT/Tables/wtpdiff_ols.tex", replace label 
		drop(0.inac_bel3* ?.inac_bel3#0* ?.risk_pref 0.risk_pref#* 1.inac_bel3#3.*#*neg *3.risk_pref*)
		order(false_pos false_neg 1.risk_pref#c.false_pos 1.risk_pref#c.false_neg 2.risk_pref#c.false_pos 2.risk_pref#c.false_neg ?.inac_bel3 ?.inac_bel3#c.* ?.inac_bel3#1.risk_pref ?.inac_bel3#1.risk_*c.*)
		cells(b(fmt(3)) se(par star fmt(3))) fragment booktabs 
		starlevels(* 0.1 ** 0.05 *** 0.01) varwidth(15) mlabels(, none) collabels(, none) nomtitle noobs nonum noline
		stat(r2 N b_RA_fp_X se_RA_fp_X b_RA_fn_X se_RA_fn_X,  label("\midrule $ R^2$" "Obs" "Risk-Averse Subjects: \\ \hspace{0.5em} False Positive" " " "\hspace{0.5em} False Negative" " ") fmt(3 0 3 3 3 3) layout(@ @ @ (@) @ (@)))
		substitute(_ \_ ) style(tex);
#delimit cr
*/

#delimit ;
	esttab * using "$ROOT/Tables/wtpdiff_ols.tex", replace label 
		drop(0.inac_bel3* ?.inac_sub 0.inac_sub#* ?.risk_pref 0.risk_pref#* *3.risk_pref*)
		order(false_pos false_neg 1.risk_pref#c.false_pos 1.risk_pref#c.false_neg 2.risk_pref#c.false_pos 2.risk_pref#c.false_neg)
		cells(b(fmt(3)) se(par star fmt(3))) fragment booktabs 
		starlevels(* 0.1 ** 0.05 *** 0.01) varwidth(15) mlabels(, none) collabels(, none) nomtitle noobs nonum noline
		stat(r2 N b_RA_fp_X se_RA_fp_X p_RA_fp_X b_RA_fn_X se_RA_fn_X p_RA_fn_X b_RL_fp_X se_RL_fp_X p_RL_fp_X b_RL_fn_X se_RL_fn_X p_RL_fn_X ,  
		label("\midrule $ R^2$" "Obs"  "[1em] Risk-Averse Subjects: \\ \hspace{0.5em} False Positive" "\hspace{1em}  se" "\hspace{1em} $ p$-value" "[0.5em] \hspace{0.5em} False Negative" "\hspace{1em}  se" "\hspace{1em}  $ p$-value" 
		"[1em] Risk-Loving Subjects: \\ \hspace{0.5em} False Positive" "\hspace{1em}  se" "\hspace{1em}  $ p$-value"  "[0.5em] \hspace{0.5em} False Negative" "\hspace{1em}  se" "\hspace{1em}  $ p$-value" ) 
		fmt(3 0 3 3 3 3 3 3 3 3 3 3 3 3) layout(@ @ @ (@) [@] @ (@) [@] @ (@) [@] @ (@) [@]))
		substitute(_ \_ ) style(tex);
#delimit cr


XXX
*eststo: reghdfe wtp_diff i.accur_bel3#i.risk_pref##c.false_pos i.accur_bel3#i.risk_pref##c.false_neg if plevel<300, abs(subject_id) vce(cluster subject_id)
esttab using "./Tables/table_wtpdiff_01rfull_ag.tex", b(%9.3g) se(%9.1f) ar2(%9.2f) label title(WTP minus Value of Information (OLS)) mtitles("" "" "" "" "") star("*" 0.10 "**" 0.05 "***" 0.01) nobaselevels compress nogaps replace
esttab using "./Tables/table_wtpdiff_01r_ag.tex", b(%9.3g) se(%9.1f) ar2(%9.2f) label title(WTP minus Value of Information (OLS)) mtitles("" "" "" "" "") star("*" 0.10 "**" 0.05 "***" 0.01) noomit nobaselevels compress nogaps replace

