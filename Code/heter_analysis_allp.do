set more off
clear all

**Requires: estout, moremata

**Doing the latent class heterogeneity analysis including all the priors

*!!put your working folder below:
cd C:\Tornado_warnings\Experiment\Alerts_Experiment


set seed 135

use "./Temp/prep_heter.dta", replace

*Labeling self-reported strategies
gen strategy_short=strategy_ip
replace strategy_short=3 if strategy_short>3
label var strategy_short "Strategy"
label define strategy_short_l 1 "Rational" 2 "Seek honest" 3 "Other"
label values strategy_short strategy_short_l



eststo clear
eststo: reg ip_ i.strategy_short##c.phintWB i.strategy_short##c.phintBW, vce(cluster subject_id)
eststo: reg ip_ i.strategy_short##c.phintWB i.strategy_short##c.phintBW if hint=="w", vce(cluster subject_id)
eststo: reg ip_ i.strategy_short##c.phintWB i.strategy_short##c.phintBW if hint=="b", vce(cluster subject_id)
esttab using "./Tables/table_ip_strat.tex", b(%9.3g) t(%9.1f) ar2(%9.2f) label addnotes("Reporting average coeffs, errors are clustered by subject.") mtitles("" "S=White" "S=Black") title(Informed protection response and reported strategy: LPM) star("*" 0.10 "**" 0.05 "***" 0.01) nobaselevels compress nogaps replace


*Identify people who follow particular strategies closely
*Hint-followers: protect when black
gen ip_black=ip_
gen ip_white=ip_
replace ip_black=. if hint=="w"
replace ip_white=. if hint=="b"
bys subject_id: egen avprot_b=mean(ip_black)
bys subject_id: egen avprot_w=mean(ip_white)

gen hint_follower=(avprot_b==1)&(avprot_w==0)
tab hint_follower

*Prior followers:
*Always protect in last 3 rounds but never in first 3
gen ip_low=ip_
gen ip_high=ip_
replace ip_low=. if round>3
replace ip_high=. if round<4
bys subject_id: egen avprot_l=mean(ip_low)
bys subject_id: egen avprot_h=mean(ip_high)
gen p_follower=(avprot_l<0.17)&(avprot_h>0.82)
tab p_follower
*There are no p followers!!

*Bayesians:
save "./Temp/strat_analysis.dta", replace
sort subject_id post_prob

*Cleaner way to find violations:
*finding non-monotonicites with respect to posterior probability in the informed protection task:
collapse (mean) ip_ avprot_b, by(subject_id post_prob)
sort subject_id post_prob
gen sort_violation=0
replace sort_violation=1 if ((ip_==0)&(ip_[_n-1]>0))&(post_prob>post_prob[_n-1])&(subject_id==subject_id[_n-1])
replace sort_violation=1 if ip_>0&post_prob==0
replace  sort_violation=1 if ip_<avprot_b&post_prob==1
bys subject_id: egen bayesian=max(sort_violation)
replace bayesian=1-bayesian
tab bayesian
drop sort_violation
collapse (mean) bayesian, by(subject_id)
save "./Temp/bayesians.dta", replace
use "./Temp/strat_analysis.dta", replace
merge m:1 subject_id using "./Temp/bayesians.dta"
drop _merge

replace bayesian=0 if ip_==0&post_prob==1&avprot_b>0
bys subject_id: egen bayesian0=min(bayesian)
replace bayesian=bayesian0
drop bayesian0
tab bayesian

*Belief-followers:
sort subject_id be_
gen sort_violation=0
replace sort_violation=1 if ((ip_==0)&(ip_[_n-1]==1))&(be_>be_[_n-1])&(subject_id==subject_id[_n-1])
replace sort_violation=1 if ip_==1&be_==0
replace  sort_violation=1 if ip_==0&be_==1&avprot_b<1
bys subject_id: egen bel_follower=max(sort_violation)
replace bel_follower=1-bel_follower
tab bel_follower
drop sort_violation


*Subjects protecting when having black-eyed on white with no white-eyed:
gen overprotect0=0
replace overprotect0=1 if phintWB==0&phintBW>0&hint=="w"&ip_==1
bys subject_id: egen overprotect=max(overprotect0)
drop overprotect0

gen underprotect0=0
replace underprotect0=1 if ip_==0&phintBW==0&phintWB>0&avprot_b>0&hint=="b"
bys subject_id: egen underprotect=max(underprotect0)
drop underprotect0

tab underprotect overprotect
tab underprotect bayesian
tab overprotect bayesian

save "./Temp/strat_analysis.dta", replace


use "./Temp/strat_analysis.dta", replace
*Just knowing the prop of dishonest gremlins is almost as good for predicting decisions as knowing FP and FN!
gen false_prob=phintWB+phintBW



*Make the dataset balanced:
drop if missing(ip_)
sort subject_id round hint
by subject_id: egen nchoices=count(round)
drop if nchoices<12


expand 2, gen(alt)
replace alt=alt+1
replace ip_=1-ip_ if alt==1 //ip_ designates choosing no protection case here//
egen choice_id=group(subject_id round hint)
sort subject_id round hint alt


order subject_id round hint choice_id alt ip_ post_prob false_prob

bro
replace alt=alt-1 //alt=1 means choosing protection now//

gen phigh=post_prob
replace phigh=0 if alt==0

gen phigh_be=be_
replace phigh_be=0 if alt==0


gen phigh_hint=blackhint
replace phigh_hint=0 if alt==0

reg phigh_hint phigh
predict phigh_hintd, residual

reg phigh_be phigh
predict phigh_bed, residual

gen phigh_hintd2=-phigh_hintd


gen phigh_acc=blackhint*(1-false_prob)
replace phigh_acc=0 if alt==0



gen phintWB1=phintWB
gen phintBW1=phintBW
*replace phintWB=0 if blackhint==0
*replace phintBW=0 if blackhint==1
replace phintWB=0 if alt==0
replace phintBW=0 if alt==0

replace false_prob=0 if alt==0
bys alt: sum phigh

gen prior=p
replace p=0 if alt==1



reg p phigh phigh_hintd false_prob
predict pdev, residual

label var p "Prior"
label var phigh "Posterior"
label var phigh_hint "Black hint"
label var phigh_hintd "Black hint (res)"
label var false_prob "Prop. of liars"
label var false_prob "Prop. of lying gremlins"
gen phintB=prior*phintBB+(1-prior)*phintBW //probabiity to observe the black signal
gen phintW=1-phintB
save "./Temp/lc_analysis.dta", replace


**Doing some simple exploratory specifications:
*lclogit2 ip_ alt, rand(p phigh_hint false_prob) group(choice_id) id(subject_id) nclasses(3) seed(10001)
*logit ip_ phigh phigh_hintd if alt==1

*logit ip_ phigh if alt==1 //matching estimates
*clogit ip_ alt phigh, group(choice_id) //matching estimates


*logit ip_ phigh i.subject_id if alt==1 //not matching estima*tes
*clogit ip_ alt phigh, group(subject_id) //not matching estimates


*logit ip_ phigh_hint phintWB phintBW if alt==1 //matching estimates 2
*clogit ip_ alt phigh_hint phintWB phintBW, group(choice_id) //matching estimates 2




***Baseline estimation***
*matrix drop class_compar
*One class: do it by imposing constraints between the two classes (as lclogit2 doesn't allow for one class)
constraint 1 [Class1]phigh_hint = [Class2]phigh_hint
constraint 2 [Class1]phigh = [Class2]phigh
constraint 3 [Class1]false_prob = [Class2]false_prob
lclogit2 ip_ alt, rand(phigh_hintd false_prob phigh) group(choice_id) id(subject_id) nclasses(2) constraints(1 2 3) seed(10001)
local akaike1=e(aic)
local nclasses=1
matrix coeff=e(b)
display(`akaike1')
local bic=e(bic)
display(`bic')
local akaike_cor=e(caic)
display(`akaike_cor')
matrix start = e(b)
matrix class_compar = nullmat(class_compar) \ `nclasses', `bic', `akaike1', `akaike_cor'
matrix colnames class_compar = "N Classes" "BIC" "AIC"  "AIC corrected"

matrix class1_coeff=coeff["y1", "Class1:"]
matrix lc_results=1, 1, coeff["y1", "Fix:alt"], class1_coeff, 1, `bic' 
matrix list lc_results

forvalues nclasses = 2/4 {
  lclogit2 ip_ alt, rand(phigh_hintd false_prob phigh) group(choice_id) id(subject_id) nclasses(`nclasses') seed(10001)
  local akaike1=e(aic)
  display(`akaike1')
  local bic=e(bic)
  display(`bic')
  local akaike_cor=e(caic)
  display(`akaike_cor')
  matrix class_compar = nullmat(class_compar) \ `nclasses', `bic', `akaike1', `akaike_cor'
  matrix start = e(b)
  lclogitml2 ip_ alt, rand(phigh_hintd false_prob phigh) group(choice_id) id(subject_id) nclasses(`nclasses')  from(start)
}
matlist class_compar, name(columns)
esttab matrix(class_compar) using "./Tables/class_search.tex", tex replace




lclogit2 ip_ alt, rand(phigh_hintd false_prob phigh) group(choice_id) id(subject_id) nclasses(2) seed(10001)
local akaike1=e(aic)
local nclasses=2
display(`akaike1')
local bic=e(bic)
display(`bic')
local akaike_cor=e(caic)
display(`akaike_cor')
matrix start = e(b)
lclogitml2 ip_ alt, rand(phigh_hintd false_prob phigh) group(choice_id) id(subject_id) nclasses(2)  from(start)
display("coefficients:")
matrix coeff=e(b)
matrix shares=e(P)
matrix list coeff

*Assembling the matrix of results:
matrix class1_coeff=coeff["y1", "Class1:"]
matrix class2_coeff=coeff["y1", "Class2:"]
local modelname=1
local classn=1
matrix lc_results = lc_results \ 2, 1, coeff["y1", "Fix:alt"], class1_coeff, shares[1,1], `bic' \ 2, 2, coeff["y1", "Fix:alt"], class2_coeff, shares[2,1], `bic'
matrix list lc_results
matrix colnames lc_results = "Model" "Class" "Alt"  "Hint" "False_prob" "Posterior" "Class share" "BIC"
matrix list lc_results


lclogit2 ip_ alt, rand(phigh_hintd false_prob phigh) group(choice_id) id(subject_id) nclasses(3) seed(10001)
local akaike1=e(aic)
local nclasses=3
display(`akaike1')
local bic=e(bic)
display(`bic')
local akaike_cor=e(caic)
display(`akaike_cor')
matrix start = e(b)
lclogitml2 ip_ alt, rand(phigh_hintd false_prob phigh) group(choice_id) id(subject_id) nclasses(3)  from(start)
display("coefficients:")
matrix coeff=e(b)
matrix shares=e(P)
matrix list coeff

*Assembling the matrix of results:
matrix class1_coeff=coeff["y1", "Class1:"]
matrix class2_coeff=coeff["y1", "Class2:"]
matrix class3_coeff=coeff["y1", "Class3:"]
local modelname=1
local classn=1
matrix lc_results = lc_results \ 3, 1, coeff["y1", "Fix:alt"], class1_coeff, shares[1,1], `bic'
matrix lc_results = lc_results  \ 3, 2, coeff["y1", "Fix:alt"], class2_coeff, shares[2,1], `bic'
matrix lc_results = lc_results  \ 3, 3, coeff["y1", "Fix:alt"], class3_coeff, shares[3,1], `bic'
matrix list lc_results
matrix colnames lc_results = "Model" "Class" "Alt"  "Hint" "False_prob" "Posterior" "Class share" "BIC"
matrix rownames lc_results=_:
matrix roweq lc_results=""
matrix rownames lc_results=""
matrix coleq lc_results = "" 
matrix list lc_results
esttab matrix(lc_results) using "./Tables/lc_results_all.tex", b(%9.3g)  sfmt(%9.3g)tex replace




*predicting posterior probabilities of belonging to each class given choices
*re-estimating for 2 classes
lclogit2 ip_ alt, rand(phigh_hintd false_prob phigh) group(choice_id) id(subject_id) nclasses(2) seed(10001)
local akaike1=e(aic)
local nclasses=2
matrix start = e(b)
lclogitml2 ip_ alt, rand(phigh_hintd false_prob phigh) group(choice_id) id(subject_id) nclasses(2)  from(start)

lclogitpr2 classpr, cp
sum classpr*

*classification confidence (max prob out of all the classes):
egen double cpmax = rowmax(classpr1-classpr2)
bys subject_id: gen first = _n==1
sum cpmax if first==1

*prior probabilities:
lclogitpr2 pr, pr

*Checking classification quality:
generate byte class = .

forvalues c = 1/`e(nclasses)' {
    quietly replace class = `c' if cpmax==classpr`c'
}

label define class_l 1 "Honesty seekers" 2 "Cautious Bayesians"
label values class class_l

tab class if first==1

forvalues c = 1/`e(nclasses)' {
   quietly summarize pr if class == `c' & ip_==1
   local n=r(N)
   local a=r(mean)
   quietly summarize classpr`c' if class == `c' & ip_==1
   local b=r(mean)
   matrix pr = nullmat(pr) \ `n', `c', `a', `b'
}

matrix colnames pr = "Obs" "Class" "Uncond_Pr" "Cond_PR"
matlist pr, name(columns)
esttab matrix(pr) using "./Tables/class_pr.tex", tex replace

tab class stat_educ if first==1, chi2

gen prot_cost=5
gen loss=20

replace false_prob=phintWB+phintBW
eststo clear
eststo: logit ip_ blackhint false_prob post_prob if alt==1&class==1, vce(cluster subject_id)
local r2p=e(r2_p)
local llike=e(ll)
predict ip_class1
eststo m1: margins, dydx(blackhint false_prob post_prob) post
estadd scalar r2p = `r2p'
estadd scalar llike = `llike'

eststo: logit ip_ blackhint false_prob post_prob if alt==1&class==2, vce(cluster subject_id)
local r2p=e(r2_p)
local llike=e(ll)
predict ip_class2
eststo m2: margins, dydx(blackhint false_prob post_prob) post
estadd scalar r2p = `r2p'
estadd scalar llike = `llike'
esttab m1 m2 using "./Tables/table_ipclasses.tex", b(%9.3g) t(%9.1f) ar2(%9.2f) stats(N r2p llike, labels("N" "Pseudo R-squared" "Log-likelihood")) label addnotes("Errors are clustered by subject, average marginal treatment effects") mtitles("Honesty Seekers" "Cautious Bayesians") title(IP response by class) star("*" 0.10 "**" 0.05 "***" 0.01) nobaselevels compress nogaps replace



keep if alt==1 //returning back to one obs per choice structure
replace p=prior
*efficiency analysis
*informed protection losses:

gen loss_class1=ip_class1*prot_cost+(1-ip_class1)*post_prob*loss
replace loss_class1=phintB*loss_class1 if blackhint==1
replace loss_class1=phintW*loss_class1 if blackhint==0


gen loss_class2=ip_class2*prot_cost+(1-ip_class2)*post_prob*loss
replace loss_class2=phintB*loss_class2 if blackhint==1
replace loss_class2=phintW*loss_class2 if blackhint==0

sum loss_class1 if class==1
sum loss_class2 if class==2
sum loss_class1 loss_class2

*calculate minimal expected losses 
gen ip_o=post_prob>(prot_cost/loss)
gen loss_opt=ip_o*prot_cost+(1-ip_o)*post_prob*loss
replace loss_opt=phintB*loss_opt if blackhint==1
replace loss_opt=phintW*loss_opt if blackhint==0
sum loss_opt

*losses relative to the switching point:
logit ip_ i.blackhint##c.phintBW i.blackhint##c.phintWB i.plevel, vce(cluster subject_id)
predict ip_base
gen loss_base=ip_base*prot_cost+(1-ip_base)*post_prob*loss
replace loss_base=phintB*loss_base if blackhint==1
replace loss_base=phintW*loss_base if blackhint==0
sum loss_base

*calculate actual losses:
gen loss_ind=ip_*prot_cost+(1-ip_)*post_prob*loss
replace loss_ind=phintB*loss_ind if blackhint==1
replace loss_ind=phintW*loss_ind if blackhint==0
sum loss_ind
bys class: sum loss_ind

save "./Temp/xtlosses.dta", replace
use "./Temp/xtlosses.dta", replace
sum ip_class1 if class==1
sum ip_class2 if class==2
collapse (sum) loss_ind loss_class1 loss_class2 loss_base loss_opt (first) class strategy_short classpr2 age sex college stat_educ accur_bel mention_prop mention_hon mention_hint strategy_ip switchprob totprot ncorrect informed_correct wtpq_correct, by(subject_id plevel phintWB phintBW)

*create the matrix of results
sum loss_class1 if class==1
local L_class1=r(mean)
display(`L_class1')
sum loss_class2 if class==2
local L_class2=r(mean)
display(`L_class2')
sum loss_base
local L_base=r(mean)
sum loss_opt
local L_opt=r(mean)
display(`L_opt')
local rat_class1=100.0*`L_class1'/`L_opt'
local rat_class2=100.0*`L_class2'/`L_opt'
local rat_base=100.0*`L_base'/`L_opt'
display(`rat_class1')
local rat_opt=100

matrix losses_matr=(`L_class1', `rat_class1' \ `L_class2', `rat_class2' \ `L_base', `rat_base' \ `L_opt', `rat_opt')
matrix colnames losses_matr = "Mean loss" "\\% of optimal"
matrix rownames losses_matr = "Honesty Seekers" "Cautious Bayesians" "Baseline (All)" "Optimal"
matrix list losses_matr
esttab matrix(losses_matr) using "./Tables/losses_matr.tex", b(%9.3g) title("Expected IP losses by strategy") sfmt(%9.3g)tex replace



*Saving the classification for later:
collapse (mean) classpr2 (first) class strategy_short age sex college stat_educ accur_bel mention_prop mention_hon mention_hint strategy_ip switchprob totprot informed_correct, by(subject_id)
save "./Temp/ipclasses.dta", replace

use "./Temp/ipclasses.dta", replace

tab class mention_prop, chi2
tab class mention_hint, chi2
tab class mention_hon, chi2


tab class strategy_short, chi2
tab class strategy_ip, chi2

tab class sex, chi2
tab class age, chi2
tab class stat_educ, chi2

reg classpr2 i.strategy_short sex stat_educ college accur_bel

label var strategy_short "Strategy"
label define strategy_short_l 1 "Rational" 2 "Seek honest" 3 "Other"
label values strategy_short strategy_short_l

label var classpr2 "Strategy2 prob"
label var switchprob "Switching prob (BP)"
label var age "Age"
label var sex "Female"
label var stat_educ "Stat. classes"
label var accur_bel "Accur. beliefs"
label var informed_correct "IP quiz"
label var totprot "RA measure0"

eststo clear
eststo: reg classpr2 i.strategy_short, vce(robust)
eststo: reg classpr2 sex age stat_educ informed_correct, vce(robust)
eststo: reg classpr2 accur_bel totprot, vce(robust)
esttab using "./Tables/table_class_det.tex", b(%9.3g) se(%9.1f) ar2(%9.2f) label title(Correlates of Strategies Used) mtitles("" "" "") star("*" 0.10 "**" 0.05 "***" 0.01) nobaselevels compress nogaps replace







*********************************************************************
**Alternative model: hint-specific sensitivity to FP and FN rates:
*drop phintWB0 phintWB1 phintBW0 phintBW1
use "./Temp/lc_analysis.dta", replace

drop phintWB1 phintBW1

gen phintWB0=phintWB
replace phintWB0=0 if phigh_hint==1

gen phintBW0=phintBW
replace phintBW0=0 if phigh_hint==1

gen phintWB1=phintWB
replace phintWB1=0 if phigh_hint==0

gen phintBW1=phintBW
replace phintBW1=0 if phigh_hint==0

label var phintWB0 "FN rate*White hint"
label var phintWB1 "FN rate*Black hint"
label var phintBW0 "FP rate*White hint"
label var phintBW1 "FP rate*Black hint"


sum phintWB0 phintWB1 phintBW0 phintBW1


sum phintWB0 phintWB1 phintBW0 phintBW1 if phigh_hint==1

gen ip_all=ip_

lclogit2 ip_ alt, rand(phigh_hint phintWB0 phintBW0 phintWB1  phintBW1) group(choice_id) id(subject_id) nclasses(2) seed(10001)
local akaike1=e(aic)
local nclasses=2
display(`akaike1')
local bic=e(bic)
display(`bic')
local akaike_cor=e(caic)
display(`akaike_cor')
matrix start = e(b)
lclogitml2 ip_ alt, rand(phigh_hint phintWB0 phintBW0 phintWB1  phintBW1) group(choice_id) id(subject_id) nclasses(2)  from(start)
display("coefficients:")
matrix coeff=e(b)
matrix shares=e(P)
matrix list coeff

gen prot_cost=5
gen loss=20
local nclasses=2

lclogitpr2 classpr, cp
sum classpr*

*classification confidence (max prob out of all the classes):
egen double cpmax = rowmax(classpr1-classpr2)
bys subject_id: gen first = _n==1
sum cpmax if first==1

*prior probabilities:
lclogitpr2 pr, pr

*Checking classification quality:
generate byte class = .

forvalues c = 1/`e(nclasses)' {
    quietly replace class = `c' if cpmax==classpr`c'
	
}

label define class_l 1 "Honesty seekers" 2 "Bayesians"
label values class class_l



forvalues c = 1/`e(nclasses)' {
   quietly summarize pr if class == `c' & ip_==1&plevel<300
   local n=r(N)
   local a=r(mean)
   quietly summarize classpr`c' if class == `c' & ip_==1&plevel<300
   local b=r(mean)
   matrix pr = nullmat(pr) \ `n', `c', `a', `b'
}




tab class stat_educ if first==1, chi2

eststo clear
eststo: logit ip_ blackhint phintWB0 phintBW0 phintWB1  phintBW1 if alt==1, vce(cluster subject_id)
*test phintBW=phintWB
*local p=r(p)
local r2p=e(r2_p)
local llike=e(ll)
predict ip_class0
eststo m1: margins, dydx(blackhint phintWB0 phintBW0 phintWB1  phintBW1) post
estadd scalar r2p = `r2p'
estadd scalar llike = `llike'

eststo: logit ip_ blackhint phintWB0 phintBW0 phintWB1  phintBW1 if alt==1&class==1, vce(cluster subject_id)
local r2p=e(r2_p)
local llike=e(ll)
predict ip_class1
eststo m2: margins, dydx(blackhint phintWB0 phintBW0 phintWB1  phintBW1) post
estadd scalar r2p = `r2p'
estadd scalar llike = `llike'

eststo: logit ip_ blackhint phintWB0 phintBW0 phintWB1  phintBW1 if alt==1&class==2, vce(cluster subject_id)
local r2p=e(r2_p)
local llike=e(ll)
predict ip_class2
eststo m3: margins, dydx(blackhint phintWB0 phintBW0 phintWB1  phintBW1) post
estadd scalar r2p = `r2p'
estadd scalar llike = `llike'
esttab m1 m2 m3 using "./Tables/table_ipclasses2.tex", b(%9.3g) t(%9.1f) ar2(%9.2f) stats(N r2p llike, labels("N" "Pseudo R-squared" "Log-likelihood")) label addnotes("Errors are clustered by subject, average marginal treatment effects") mtitles("All" "Class 1" "Class 2") title(IP response by class) star("*" 0.10 "**" 0.05 "***" 0.01) nobaselevels compress nogaps replace




keep if alt==1 //returning back to one obs per choice structure
*efficiency analysis
*informed protection losses:

gen loss_class0=ip_class0*prot_cost+(1-ip_class0)*post_prob*loss
replace loss_class0=phintB*loss_class0 if blackhint==1
replace loss_class0=phintW*loss_class0 if blackhint==0

gen loss_class1=ip_class1*prot_cost+(1-ip_class1)*post_prob*loss
replace loss_class1=phintB*loss_class1 if blackhint==1
replace loss_class1=phintW*loss_class1 if blackhint==0


*Also calculate N of violations of EU when protecting in response to post_prob1, but not protecting in response to post_prob2>post_prob1:

gen loss_class2=ip_class2*prot_cost+(1-ip_class2)*post_prob*loss
replace loss_class2=phintB*loss_class2 if blackhint==1
replace loss_class2=phintW*loss_class2 if blackhint==0

sum loss_class1 if class==1
sum loss_class2 if class==2
sum loss_class1 loss_class2

*calculate minimal expected losses 
gen ip_o=post_prob>(prot_cost/loss)
gen loss_opt=ip_o*prot_cost+(1-ip_o)*post_prob*loss
replace loss_opt=phintB*loss_opt if blackhint==1
replace loss_opt=phintW*loss_opt if blackhint==0
sum loss_opt

*calculate actual losses:
gen loss_ind=ip_all*prot_cost+(1-ip_all)*post_prob*loss
replace loss_ind=phintB*loss_ind if blackhint==1
replace loss_ind=phintW*loss_ind if blackhint==0
sum loss_ind
bys class: sum loss_ind

save "./Temp/xtlosses2.dta", replace
use "./Temp/xtlosses2.dta", replace
sum ip_class1 if class==1
sum ip_class2 if class==2

*calculate risk of severe loss by class:
gen Lrisk_class1=(1-ip_class1)*post_prob
replace Lrisk_class1=phintB*Lrisk_class1 if blackhint==1
replace Lrisk_class1=phintW*Lrisk_class1 if blackhint==0

gen Lrisk_class2=(1-ip_class2)*post_prob
replace Lrisk_class2=phintB*Lrisk_class2 if blackhint==1
replace Lrisk_class2=phintW*Lrisk_class2 if blackhint==0

gen Lrisk_opt=(1-ip_o)*post_prob
replace Lrisk_opt=phintB*Lrisk_opt if blackhint==1
replace Lrisk_opt=phintW*Lrisk_opt if blackhint==0

gen Lrisk0=(1-ip_class0)*post_prob
replace Lrisk0=phintB*Lrisk0 if blackhint==1
replace Lrisk0=phintW*Lrisk0 if blackhint==0

sum Lrisk0 Lrisk_class1 Lrisk_class2 Lrisk_opt
sum Lrisk_class1 if class==1
sum Lrisk_class2 if class==2

drop highprob
gen highprob=plevel>200
tab highprob

save "./Temp/loss_analysis.dta", replace

collapse (mean) loss_class1 loss_class2 loss_class0 loss_opt Lrisk0 Lrisk_class1 Lrisk_class2 Lrisk_opt, by(highprob)
bro
rename Lrisk0 Lrisk_class0
local losses loss_class0 loss_class1 loss_class2

forvalues  i=0/2 {
 display(`i')
 gen rat_loss`i'=100*loss_class`i'/loss_opt
 *gen rat_Lrisk`i'=100*Lrisk_class`i'/Lrisk_opt
}

bro
order highprob loss_class0 loss_class1 loss_class2 loss_opt rat_loss0 rat_loss1 rat_loss2 Lrisk_class0 Lrisk_class1 Lrisk_class2 Lrisk_opt
**Convert data into a matrix
mkmat *, mat(losses_matr)
matrix losses_matr=losses_matr'
matrix panel1=(losses_matr[2..5,1], (losses_matr[6..8,1]\1), losses_matr[9..12,1])
matrix panel2=(losses_matr[2..5,2], (losses_matr[6..8,2]\1), losses_matr[9..12,2])
matrix losses_matr=(panel1, panel2)
matrix list losses_matr



matrix colnames losses_matr = "Mean loss" "\\% of optimal" "Loss prob." "Mean loss" "\\% of optimal" "Loss prob."
matrix rownames losses_matr = "Baseline (all)" "Honesty seekers" "Bayesians"  "Optimal"
matrix list losses_matr
esttab matrix(losses_matr) using "./Tables/losses_matr2.tex", b(2) t(2) p(2) r2(2) title("Expected IP losses by strategy") sfmt(2) tex replace

*esttab matrix(losses_matr) using "./Tables/losses_matr2.tex", cell(fmt(2)) title("Expected IP losses by strategy") tex replace
use "./Temp/loss_analysis.dta", replace

gen fp_env=phintBW>0
gen fn_env=phintWB>0

collapse (mean) loss_ind loss_class1 loss_class2 loss_class0 loss_opt Lrisk0 Lrisk_class1 Lrisk_class2 Lrisk_opt, by(plevel fp_env fn_env)
gen los_rat=loss_ind/loss_opt
bro


use "./Temp/loss_analysis.dta", replace
*Saving the classification for later:
collapse (mean) classpr2 (first) class strategy_short age sex college stat_educ accur_bel mention_prop mention_hon mention_hint strategy_ip switchprob totprot ncorrect informed_correct wtpq_correct, by(subject_id)
save "./Temp/ipclasses2.dta", replace
use "./Temp/ipclasses2.dta", replace


tab class mention_prop, chi2
tab class mention_hint, chi2
tab class mention_hon, chi2


tab class strategy_short, chi2
tab class strategy_ip, chi2

tab class sex, chi2
tab class age, chi2
tab class stat_educ, chi2

reg classpr2 i.strategy_short sex stat_educ college accur_bel

label var strategy_short "Strategy"
label define strategy_short_l 1 "Rational" 2 "Seek honest" 3 "Other"
label values strategy_short strategy_short_l

label var classpr2 "Strategy2 prob"
label var switchprob "Switching prob (BP)"
label var age "Age"
label var sex "Female"
label var stat_educ "Stat. classes"
label var accur_bel "Accur. beliefs"
label var totprot "RA measure0"
label var informed_correct "IP quiz"

eststo clear
eststo: reg classpr2 i.strategy_short, vce(robust)
eststo: reg classpr2 sex age stat_educ, vce(robust)
eststo: reg classpr2 accur_bel totprot informed_correct, vce(robust)
esttab using "./Tables/table_class_det2.tex", b(%9.3g) se(%9.1f) ar2(%9.2f) label title(Correlates of Strategies Used) mtitles("" "" "") star("*" 0.10 "**" 0.05 "***" 0.01) nobaselevels compress nogaps replace

*First classification vs second classification:
rename class class2
rename classpr2 classpr2new
merge 1:1 subject_id using "./Temp/ipclasses.dta"
tab class class2

