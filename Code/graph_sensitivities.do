**ANALYZE THE FIRST WAVE (planned 100 participants)***
set more off
clear all

cd C:\Tornado_warnings\Experiment\Alerts_Experiment
use "./Temp/wtp_discrepancy0.dta", replace
drop highprob
gen highprob=p>0.201

sum value wtp


**getting and storing WTP regression estimates (sensitivities)
tempname p1
postfile `p1' wtp_FP wtp_FN plevel using "./Temp/sensit_estimates.dta", replace


foreach plev of numlist 100 200 300 500{
    *tobit wtp phintBW phintWB if plevel==`plev', ll(0) ul(5)
	reg wtp phintBW phintWB if plevel==`plev'
	local wtp_FP=-_b[phintBW]
	local wtp_FN=-_b[phintWB]
	post `p1' (`wtp_FP') (`wtp_FN') (`plev')
  }

postclose `p1'



**creating a frame of p values
clear all
set obs 54
gen ind=_n
gen ind_p=floor((ind-1)/9)
gen ind_fp=floor((mod(ind-1,9))/3)
gen ind_fn=mod(ind-1,3)
gen p=0.1*ind_p
recode ind_fp (0 = 0) (1 = 0.2) (2 = 0.333333333), gen(phintBW)
recode ind_fn (0 = 0) (1 = 0.2) (2 = 0.333333333), gen(phintWB)
drop if p==0


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

*dropping combinations not present in the experimental data:
drop if ind_fp+ind_fn==3
bys p: sum value

gen plevel=round(1000*p)

*adding observations to make it comparable with the WTP analysis
expand 3 if phintWB==0&phintBW==0

**getting and storing WTP regression estimates (sensitivities)
tempname p1
postfile `p1' val_FP val_FN plevel using "./Temp/sensit_estimates_val.dta", replace


*Note: ML finds a flat region (not enough degrees of freedom?)
foreach plev of numlist 100 200 300 400 500{
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

sum sensFP sensFN
label var sensFP "FP sensitivity (theory)"
label var sensFN "FN sensitivity (theory)"
label var val_FP "FP sensitivity (theory)"
label var val_FN "FN sensitivity (theory)"
label var wtp_FP "FP sensitivity (empirics)"
label var wtp_FN "FN sensitivity (empirics)"
* title("Theoretical and Empirical WTP" "Sensitivities to FP and FN rates")

*net install scheme-modern, from("https://raw.githubusercontent.com/mdroste/stata-scheme-modern/master/")
set scheme modern
#delimit ;
twoway (line val_FP p, lcolor(dknavy) lwidth(thick)) (line val_FN p, lcolor(red) lwidth(thick))
 (scatter wtp_FP p, mcolor(dknavy)) (scatter wtp_FN p, mcolor(red)),
 xtitle("Prior probability")
 legend(ring(0) position(10) bmargin(large) cols(2))
 note("OLS estimates of sensitivity to FP and FN rates by prior probability of a black ball.");
 
#delimit cr
graph export "./Graphs/sensit_comparison.png", width(1000) height(600) replace
 

