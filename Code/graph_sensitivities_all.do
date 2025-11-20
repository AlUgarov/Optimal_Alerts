**ANALYZE THE FIRST WAVE (planned 100 participants)***
set more off
clear all

cd C:\Tornado_warnings\Experiment\Alerts_Experiment
use "./Temp/wtp_discrepancy0.dta", replace
drop highprob
gen highprob=p>0.201

sum value wtp
*collapse (mean) wtp, by(phintBW phintWB p)
*stop

*stop

**getting and storing WTP regression estimates (sensitivities)
tempname p1
postfile `p1' wtp_FP wtp_FN wtp_FP_se wtp_FN_se plevel using "./Temp/sensit_estimates.dta", replace


foreach plev of numlist 100 200 300 500{
    *tobit wtp phintBW phintWB if plevel==`plev', ll(0) ul(5)
	reg wtp phintBW phintWB if plevel==`plev'
	local wtp_FP=-_b[phintBW]
	local wtp_FP_se=_se[phintBW]
	local wtp_FN=-_b[phintWB]
	local wtp_FN_se=_se[phintWB]
	post `p1' (`wtp_FP') (`wtp_FN') (`wtp_FP_se') (`wtp_FN_se') (`plev')
  }

postclose `p1'


tempname p1
postfile `p1' wtp_FPt wtp_FNt wtp_FP_set wtp_FN_set plevel using "./Temp/sensit_estimates_tobit.dta", replace


foreach plev of numlist 100 200 300 500{
    tobit wtp phintBW phintWB if plevel==`plev', ll(0) ul(5)
	local wtp_FPt=-_b[phintBW]
	local wtp_FP_set=_se[phintBW]
	local wtp_FNt=-_b[phintWB]
	local wtp_FN_set=_se[phintWB]
	post `p1' (`wtp_FPt') (`wtp_FNt') (`wtp_FP_set') (`wtp_FN_set') (`plev')
  }

postclose `p1'







collapse (mean) wtp (count) nobs=wtp, by(p phintWB phintBW)


gen phintBB=1-phintWB
gen phintWW=1-phintBW
gen protectioncost=5
gen loss=20


gen cost_bp=min(p*loss, protectioncost) //expected blind protection cost
gen false_pos=(1-p)*phintBW*protectioncost //expected protection costs for false positives
gen true_pos=p*phintBB*protectioncost //expected protection costs for true positives

gen false_neg=p*phintWB*loss //expected false positive costs
gen cost_ip=false_neg+false_pos+true_pos //total informed protection costs
gen pos_cost=false_pos+true_pos //total protection costs in response to positive signals
gen value=max(0, cost_bp-cost_ip) //theoretical value for a risk-neutral subject

gen sensFP=(1-p)*protectioncost
gen sensFN=p*(loss-protectioncost)

replace sensFP=0 if value==0
replace sensFN=0 if value==0
hist value

bys p: sum value

gen plevel=round(1000*p)

*adding observations to make it comparable with the WTP analysis
expand 3 if phintWB==0&phintBW==0

**getting and storing WTP regression estimates (sensitivities)
tempname p1
postfile `p1' val_FP val_FN plevel using "./Temp/sensit_estimates_val.dta", replace


*Note: ML finds a flat region (not enough degrees of freedom?)
foreach plev of numlist 100 200 300 500{
    *tobit value phintBW phintWB if plevel==`plev', ll(0) ul(5)
	reg value phintBW phintWB if plevel==`plev'
	local val_FP=-_b[phintBW]
	local val_FN=-_b[phintWB]
	post `p1' (`val_FP') (`val_FN') (`plev')
  }

postclose `p1'







collapse (mean) sensFP sensFN value (first) plevel, by(p)

*use "./Temp/sensit_estimates_val.dta", replace

merge 1:1 plevel using "./Temp/sensit_estimates.dta"
drop _merge
merge 1:1 plevel using "./Temp/sensit_estimates_val.dta"
drop _merge
merge 1:1 plevel using "./Temp/sensit_estimates_tobit.dta"


sum sensFP sensFN
label var sensFP "FP sensitivity (theory)"
label var sensFN "FN sensitivity (theory)"
label var val_FP "FP sensitivity (theory)"
label var val_FN "FN sensitivity (theory)"
label var wtp_FP "FP sensitivity (empirics)"
label var wtp_FN "FN sensitivity (empirics)"
label var wtp_FPt "FP sensitivity (empirics)"
label var wtp_FNt "FN sensitivity (empirics)"
* title("Theoretical and Empirical WTP" "Sensitivities to FP and FN rates")

gen wtp_FP_lb=wtp_FP-1.96*wtp_FP_se
gen wtp_FP_ub=wtp_FP+1.96*wtp_FP_se
gen wtp_FN_lb=wtp_FN-1.96*wtp_FN_se
gen wtp_FN_ub=wtp_FN+1.96*wtp_FN_se


gen p_FP = p - 0.005
gen p_FN = p + 0.005

*net install scheme-modern, from("https://raw.githubusercontent.com/mdroste/stata-scheme-modern/master/")
set scheme modern
#delimit ;
twoway (line val_FP p, lcolor(dknavy) lwidth(thick)  lpattern(".-")) (line val_FN p, lcolor(red) lwidth(thick) lpattern(dash))
 (scatter wtp_FP p, mcolor(dknavy)  msize(medlarge)) (scatter wtp_FN p, mcolor(red) msymbol(d) msize(medlarge)) (rcap wtp_FP_lb wtp_FP_ub p_FP, lcolor(dknavy) msize(medlarge)) (rcap wtp_FN_lb wtp_FN_ub p_FN, lcolor(red) msize(medlarge)),
 xtitle("Prior probability")
 yscale(r(0 9))
 xlabel(#5, labsize(medium))
 ylabel(#5, labsize(medium))
 legend(ring(0) position(11) bmargin(large) cols(2) size(medium))
 note("OLS estimates of sensitivity to FP and FN rates by prior probability of a black ball.", size(medium));
 
#delimit cr
graph export "./Graphs/sensit_comparison.png", width(1000) height(600) replace
 
drop wtp_FP_lb wtp_FP_ub wtp_FN_lb wtp_FN_ub
gen wtp_FP_lb=wtp_FPt-1.96*wtp_FP_set
gen wtp_FP_ub=wtp_FPt+1.96*wtp_FP_set
gen wtp_FN_lb=wtp_FNt-1.96*wtp_FN_set
gen wtp_FN_ub=wtp_FNt+1.96*wtp_FN_set

set scheme modern
#delimit ;
twoway (line sensFP p, lcolor(dknavy) lwidth(thick)  lpattern(".-")) (line sensFN p, lcolor(red) lwidth(thick) lpattern(dash))
 (scatter wtp_FPt p, mcolor(dknavy)  msize(medlarge)) (scatter wtp_FNt p, mcolor(red) msymbol(d) msize(medlarge)) (rcap wtp_FP_lb wtp_FP_ub p_FP, lcolor(dknavy) msize(medlarge)) (rcap wtp_FN_lb wtp_FN_ub p_FN, lcolor(red) msize(medlarge)),
 xtitle("Prior probability")
 yscale(r(0 9))
 xlabel(#5, labsize(medium))
 ylabel(#5, labsize(medium))
 legend(ring(0) position(11) bmargin(large) cols(2) size(medium))
 note("Estimates of sensitivity to FP and FN rates by prior probability of a black ball (tobit)", size(medium));
 
#delimit cr
graph export "./Graphs/sensit_comparison_tobit.png", width(1000) height(600) replace

