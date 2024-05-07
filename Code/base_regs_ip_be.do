
**********************************************
****----MAIN REGRESSIONS (beliefs and informed protection)---***********
**********************************************
set more off
clear all
graph drop _all

*!!put your working folder below:
cd C:\Tornado_warnings\Experiment\Alerts_Experiment
*cd C:\Tornado_warnings\Optimal_Alerts

use "./Temp/base_main_waves.dta", replace
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

lpoly ip_ post_prob, bwidth(0.1) ci noscatter title("Informed Protection Response") lineopt(lwidth(1)) mlwidth(thick) xtitle("Posterior probability of a black ball") ytitle("Proportion of protection choices") legend(ring(0) bplacement(seast)) ysc(r(0 1)) ylabel(0(0.2)1)
graph export "./Graphs/ip_response_lpoly.png", width(1000) height(1000) replace

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
list subject_id totdisagr1 totdisagr2 ip_dev tot_diff if ip_==0&post_prob==1

*Notes: nonsensical responses mostly come from subjects choosing different answers to different probabilities, 
* most subjects choosing no protection when prob=1 do not choose protection when prob=0
* seems to be that these subjects are just more random: sum of deviations btw posteriors and responses is much higher for them than for others even excluding that one mistake

collapse (mean) mean=ip_ (count) n=xid, by(post_probi)
gen se=sqrt(mean*(1-mean)/n)

gen task_type="Informed"

gen post_prob=0.001*post_probi
append using "./Temp/bp_graph_data.dta"
serrbar mean se post_prob if task_type=="Informed", scale (1.96) title("Informed Protection Response") lwidth(thick) xtitle("Posterior probability of a black ball") saving(inform_prot, replace)
graph export "./Graphs/ip_response.png", width(1000) height(1000) replace


gen low=mean-1.96*se
gen high=mean+1.96*se
twoway (rcap high low post_prob if task_type=="Blind", lwidth(thick)) (rcap high low post_prob if task_type=="Informed"), title("Protection Response") xtitle("Posterior probability of a black ball") ytitle("Proportion of protection choices") legend(label(1 "Blind") label(2 "Informed"))
graph export "./Graphs/ip_response_comp.png", width(1200) height(800) replace

serrbar mean se post_prob if task_type=="Blind", scale (1.96) title("Blind Protection Response") lwidth(thick) xtitle("Probability of a black ball") ytitle("Proportion of protection choices") ysc(r(0 1)) ylabel(0(0.2)1.0) saving(blindprot, replace)

graph combine blindprot.gph inform_prot.gph, ycommon

replace high=1 if high>1
replace low=0 if low<0
twoway (rcap high low post_prob if task_type=="Informed", mlwidth(thick)), title("Protection Response") xtitle("Posterior probability of a black ball") ytitle("Proportion of protection choices")  note("The bars show 95% confidence intervals for the mean proportion of subjects " "choosing protection at each probability.")
graph export "./Graphs/ip_response.png", width(1000) height(1000) replace



use "./Temp/long_ip_dat.dta", replace
set more off


****-- INFORMED PROTECTION --****
*expected response by prior prob:
gen highprob=p>0.201
label var highprob "p$>$0.2"
label define highprob_l 0 "p$\leq$0.2" 1 "p$>$0.2"
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

label var high_phintBW "FP rate x (p$>$0.2)"
label var high_phintWB "FN rate x (p$>$0.2)"
label var black_phintBW "FP rate x (S=Black)"
label var black_phintWB "FN rate x (S=Black)"

label var white_phintBW "FP rate x (S=White)"
label var white_phintWB "FN rate x (S=White)"

label var FPFN "FP rate x FN rate"


xtset subject_id
label var phintBW "FP rate"
label var phintWB "FN rate"


**Baseline IP table:


eststo clear
eststo: reg ip_ phintBW phintWB i.subject_id highprob, vce(cluster subject_id)
test phintBW=phintWB
estadd scalar p = r(p)
estadd scalar r2a = e(a_r)
eststo: reg ip_ phintBW phintWB  highprob i.subject_id if blackhint==0, vce(cluster subject_id)
test phintBW=phintWB
estadd scalar p = r(p)
estadd scalar r2a = e(a_r)
eststo: reg ip_ phintBW phintWB highprob i.subject_id if blackhint==1, vce(cluster subject_id)
test phintBW=phintWB
estadd scalar p = r(p)
estadd scalar r2a = e(a_r)
eststo: reg ip_ phintBW phintWB highprob high_phintBW high_phintWB i.subject_id, vce(cluster subject_id)
test (phintBW=phintWB) (high_phintBW=high_phintWB)
estadd scalar p = r(p)
estadd scalar r2a = e(a_r)
eststo: reg ip_ phintBW phintWB highprob high_phintBW high_phintWB i.subject_id if blackhint==0, vce(cluster subject_id)
test (phintBW=phintWB) (high_phintBW=high_phintWB)
estadd scalar p = r(p)
estadd scalar r2a = e(a_r)
eststo: reg ip_ phintBW phintWB highprob high_phintBW high_phintWB i.subject_id if blackhint==1, vce(cluster subject_id)
test (phintBW=phintWB) (high_phintBW=high_phintWB)
estadd scalar p = r(p)
estadd scalar r2a = e(a_r)
esttab using "./Tables/table_ip0_lin.tex", b(%9.3g) t(%9.1f) ar2(%9.2f) label stats(p N r2a, labels("P(FP rate $\neq$ FN rate)" "Observations" "Adjusted \(R^{2}\)")) addnotes("Errors are clustered by subject") indicate(Subject FE = *.subject_id) mtitles("All" "S=White" "S=Black" "All" "S=White" "W=Black") title(Informed protection response: linear probability regression) star("*" 0.10 "**" 0.05 "***" 0.01) nobaselevels compress nogaps replace



*IP anomalies basic: with flexible control of posteriors:
eststo clear
eststo: logit ip_ post_prob0* phintBW phintWB highprob blackhint, vce(cluster subject_id)
eststo m1: margins, dydx(phintBW phintWB highprob blackhint) post
eststo: logit ip_ post_prob0* phintBW phintWB highprob blackhint black_phintBW black_phintWB i.subject_id, vce(cluster subject_id)
eststo m2: margins, dydx(phintBW phintWB highprob blackhint black_phintBW black_phintWB) post

eststo: logit ip_ post_prob0* phintBW phintWB  highprob blackhint black_phintBW black_phintWB i.subject_id, vce(cluster subject_id)
eststo m3: margins, dydx(phintBW phintWB highprob blackhint black_phintBW black_phintWB) post

eststo: logit ip_ post_prob0* phintBW phintWB  highprob blackhint black_phintBW black_phintWB i.subject_id, vce(cluster subject_id)
eststo m3: margins, dydx(phintBW phintWB highprob blackhint black_phintBW black_phintWB) post

eststo: logit ip_ post_prob0* phintBW phintWB highprob blackhint high_phintBW high_phintWB i.subject_id, vce(cluster subject_id)
eststo m3: margins, dydx(phintBW phintWB highprob blackhint high_phintBW high_phintWB)  post

eststo: logit ip_ post_prob0* phintBW phintWB highprob blackhint black_phintBW black_phintWB  high_phintBW high_phintWB i.subject_id, vce(cluster subject_id)
eststo m4: margins, dydx(phintBW phintWB  highprob blackhint black_phintBW black_phintWB high_phintBW high_phintWB) post
esttab m1 m2 m3 m4 using "./Tables/table_ip5.tex", b(%9.3f) t(%9.1f) ar2(%9.2f) label addnotes("Reporting average marginal effects, subject FE, errors are clustered by subject." "With flexible controls of posterior probability" ///
 ) mtitles("" "" "" "") title(Informed Protection Response: logit with flexible control for posteriors) star("*" 0.10 "**" 0.05 "***" 0.01) nobaselevels compress nogaps replace

 
*Saving the file to do heterogeneity analysis with respect to IP strategies:
save "./Temp/prep_heter.dta", replace


use "./Temp/prep_heter.dta", replace
*Test average response vs posterior probability by 



*Is it just an error in beliefs?
*Doing the IP regression and controlling for beliefs:
gen bes=be_^2
mkspline bespline1 0.2 bespline2 0.4 bespline3 0.6 bespline4 0.8 bespline5 = be_

  
set more off
**Robustness: flexible control both for beliefs and posteriors:
eststo clear:
eststo: logit ip_ post_prob0* white_phintBW white_phintWB  highprob blackhint black_phintBW black_phintWB i.subject_id, vce(cluster subject_id)
local r2p=e(r2_p)
local llike=e(ll)
eststo m1: margins, dydx(blackhint white_phintBW white_phintWB black_phintBW black_phintWB highprob) post
estadd scalar r2p = `r2p'
estadd scalar llike = `llike'

eststo: logit ip_ post_prob0* white_phintBW white_phintWB highprob blackhint black_phintBW black_phintWB  high_phintBW high_phintWB i.subject_id, vce(cluster subject_id)
local r2p=e(r2_p)
local llike=e(ll)
eststo m2: margins, dydx(blackhint white_phintBW white_phintWB  black_phintBW black_phintWB highprob high_phintBW high_phintWB) post
estadd scalar r2p = `r2p'
estadd scalar llike = `llike'


eststo: logit ip_ post_prob0* bespline* white_phintBW white_phintWB highprob blackhint black_phintBW black_phintWB i.subject_id, vce(cluster subject_id)
local r2p=e(r2_p)
local llike=e(ll)
eststo m3: margins, dydx(blackhint white_phintBW white_phintWB black_phintBW black_phintWB highprob) post
estadd scalar r2p = `r2p'
estadd scalar llike = `llike'

eststo: logit ip_ post_prob0* bespline* white_phintBW white_phintWB highprob blackhint black_phintBW black_phintWB high_phintBW high_phintWB i.subject_id, vce(cluster subject_id)
local r2p=e(r2_p)
local llike=e(ll)
eststo m4: margins, dydx(blackhint white_phintBW white_phintWB  black_phintBW black_phintWB highprob high_phintBW high_phintWB) post
estadd scalar r2p = `r2p'
estadd scalar llike = `llike'

esttab m1 m2 m3 m4 using "./Tables/table_ip_flex.tex", b(%9.3f) t(%9.3f) label addnotes("With flexible controls of posterior probability and beliefs" ///
  "Subject FE, errors are clustered by subject, average marginal treatment effects") stats(N r2p llike, labels("N" "Pseudo R-squared" "Log-likelihood") fmt(0 3 3)) mtitles("Posterior only" "Posterior only" "Both" "Both") title(Informed Protection Response: flexible control for posteriors and beliefs) star("*" 0.10 "**" 0.05 "***" 0.01) nobaselevels compress nogaps replace


  
 **Robustness: LINEAR! flexible control both for beliefs and posteriors:
eststo clear
eststo: reg ip_ post_prob0* white_phintBW white_phintWB black_phintBW black_phintWB blackhint highprob i.subject_id, vce(cluster subject_id)
eststo: reg ip_ post_prob0* white_phintBW white_phintWB black_phintBW black_phintWB highprob blackhint high_phintBW high_phintWB i.subject_id, vce(cluster subject_id)
eststo: reg ip_ post_prob0* bespline* white_phintBW white_phintWB black_phintBW black_phintWB highprob blackhint i.subject_id, vce(cluster subject_id)
eststo: reg ip_ post_prob0* bespline* white_phintBW white_phintWB black_phintBW black_phintWB highprob blackhint high_phintBW high_phintWB i.subject_id, vce(cluster subject_id)
esttab using "./Tables/table_ip_flexlin.tex", b(%9.3g) t(%9.1f)  ar2(%9.2f) label addnotes("With flexible controls of posterior probability and beliefs" ///
  "Subject FE, errors are clustered by subject") drop(_cons) indicate("Subject FE = *.subject_id" "Posterior=post_prob0*" "Beliefs=bespline*") mtitles("" "" "" "") title(Informed Protection Response: flexible control for posteriors and beliefs, LPM) star("*" 0.10 "**" 0.05 "***" 0.01) nobaselevels compress nogaps replace
esttab using "./Tables/table_ip_flexlin_pres.tex", b(%9.3g) t(%9.1f)  ar2(%9.2f) label drop(_cons) indicate("Subject FE = *.subject_id" "Posterior=post_prob0*" "Beliefs=bespline*") mtitles("" "" "" "") star("*" 0.10 "**" 0.05 "***" 0.01) nobaselevels compress nogaps replace


  


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

*Final robustness check: semiparametric control for posteriors
eststo clear
eststo: semipar ip_ phintBW phintWB, nonpar(post_prob)
eststo: semipar ip_ highprob phintBW highprobBW phintWB highprobWB, nonpar(post_prob) 
eststo: semipar ip_ blackhint phintBW phintWB blackhintBW blackhintWB, nonpar(post_prob)
eststo: semipar ip_ stat_educ phintBW phintWB statBW statWB, nonpar(post_prob)
esttab using "./Tables/table_ip5_semi.tex", b(%9.3g) t(%9.1f) ar2(%9.2f) label mtitles("" "" "" "") title(Informed protection response: semiparametric control for posteriors) star("*" 0.10 "**" 0.05 "***" 0.01) nobaselevels compress nogaps replace

replace bel_err=-bel_err

save "./Temp/prep_beliefs.dta", replace

use "./Temp/prep_beliefs.dta", replace
**Comparison of protection by information 
gen fp_env=phintBW>0
gen fn_env=phintWB>0
bys fp_env fn_env blackhint: sum ip_ 
gen prot_cost=5
gen loss=20
gen ip_o=post_prob>(prot_cost/loss)

tempname p1
postfile `p1' nrow signal false_pos false_neg prot posterior mean_opt ptest2 using "./Temp/ip_by_environment.dta", replace

forvalues k=0/1{
  forvalues i=0/1{
  
    forvalues j=0/1{
	  if `k'==0{
	    ttest ip_==0 if fp_env==`i'&fn_env==`j'&blackhint==0, level(95)
		local ptest=r(p_u)
		
	  }
	  else {
	    ttest ip_==1 if fp_env==`i'&fn_env==`j'&blackhint==1, level(95)
	  	local ptest=r(p_l)
	  }
	  local prot=r(mu_1)
	  sum post_prob if fp_env==`i'&fn_env==`j'&blackhint==`k'
	  local posterior=r(mean)
	  ttest ip_==ip_o if fp_env==`i'&fn_env==`j'&blackhint==`k'
	  local mean_opt=r(mu_2)
	  local ptest2=r(p)
	  local nrow=4*`k'+2*`i'+`j'+1
	  post `p1' (`nrow') (`k') (`i') (`j') (`prot') (`posterior') (`mean_opt') (`ptest2')
    }
  }
}

postclose `p1'

gen testgroup=0 if fp_env==0&fn_env==0&blackhint==0
replace testgroup=1 if fp_env==1&fn_env==0&blackhint==0

ttest ip_, by(testgroup)


use "./Temp/ip_by_environment.dta", replace
tostring nrow, replace
replace nrow="("+nrow+")"

tostring false_pos false_neg signal, replace
foreach var of varlist false_pos false_neg{
  replace `var'="No" if `var'=="0"
  replace `var'="Yes" if `var'=="1"
}
replace signal="White" if signal=="0"
replace signal="Black" if signal=="1"
bro
format prot posterior mean_opt ptest2 %9.3f
listtex using "./Tables/bigpicture_IP.tex", type rstyle(tabular) head("\begin{table}[H]\centering \footnotesize \caption{Average Protection by Signal Type} \begin{tabular}{cccccccc} \hline \hline" `"\textbf{Row}&\textbf{Signal}&\textbf{False-pos.}&\textbf{False-neg.}&\textbf{\% protect}& \textbf{Posterior} & \textbf{Optimal} & \textbf{P(=optimal)} \\ \hline"') foot("\hline \end{tabular} \end{table}") replace


order nrow signal false_pos false_neg posterior prot mean_opt ptest2
listtex using "./Tables/bigpicture_IP_AU.tex", type rstyle(tabular) replace

drop nrow
listtex using "./Tables/bigpicture_IP_pres.tex", type rstyle(tabular) head("\begin{table}[H]\centering \scriptsize \begin{tabular}{ccccccc} \hline \hline" `"\textbf{Signal}&\textbf{False-pos.}&\textbf{False-neg.}&\textbf{\% protect}& \textbf{Posterior} & \textbf{Optimal} & \textbf{P(=optimal)} \\ \hline"') foot("\hline \end{tabular} \end{table}") replace



****-- BELIEF ELICITATION --****

*Beliefs accuracy by signal characteristics:
*bel_err=post_prob-be_

*save "./Temp/prep_beliefs.dta", replace
use "./Temp/prep_beliefs.dta", replace
gen fp_env=phintBW>0
gen fn_env=phintWB>0

collapse (mean) ip_ post_prob, by(fp_env fn_env blackhint)
sort fp_env fn_env blackhint

*Anchoring effects for beliefs
*checking if a subject reports the same beliefs for first and second prior:
use "./Temp/prep_beliefs.dta", replace
sort subject_id round blackhint

gen be_p_change=be_~=be_[_n-6]
replace be_p_change=0 if round<4
by subject_id: egen nchanges_p=sum(be_p_change)

gen be_s_change=be_~=be_[_n-2]
replace be_s_change=0 if inlist(round,1,2,4)
by subject_id: egen nchanges_s=sum(be_s_change)
tab nchanges_s nchanges_p



use "./Temp/prep_beliefs.dta", replace
gen fp_env=phintBW>0
gen fn_env=phintWB>0




tempname p1
postfile `p1' nrow false_pos false_neg signal posterior bel_err ptest using "./Temp/bel_by_environment.dta", replace

forvalues k=0/1{
  forvalues i=0/1{
  
    forvalues j=0/1{
		sum post_prob if fp_env==`i'&fn_env==`j'&blackhint==`k'
	    local posterior=r(mean)
	    ttest bel_err==0 if fp_env==`i'&fn_env==`j'&blackhint==`k', level(95)
		local ptest=r(p)
		local bel_err=r(mu_1)
		local nrow=4*`k'+2*`i'+`j'+1
		post `p1' (`nrow') (`i') (`j') (`k') (`posterior') (`bel_err') (`ptest')
	  }
    }
  }

postclose `p1'

use "./Temp/bel_by_environment.dta", replace
tostring nrow, replace
replace nrow="("+nrow+")"

tostring false_pos false_neg signal, replace
foreach var of varlist false_pos false_neg{
  replace `var'="No" if `var'=="0"
  replace `var'="Yes" if `var'=="1"
}
replace signal="White" if signal=="0"
replace signal="Black" if signal=="1"
bro
format posterior bel_err ptest %9.3f
listtex using "./Tables/bigpicture_bel.tex", type rstyle(tabular) head("\begin{table}[H]\centering \caption{Average Belief Error by Signal Type} \begin{tabular}{ccccccc} \hline \hline" `" &\textbf{False-pos.}&\textbf{False-neg.}&\textbf{Signal}&\textbf{Posterior}&\textbf{Belief error}& \textbf{P($=0$)}\\ \hline"') foot("\hline \end{tabular} \end{table}") replace

order nrow false_pos false_neg signal posterior bel_err ptest
listtex using "./Tables/bigpicture_bel_AU.tex", type rstyle(tabular) replace



use "./Temp/prep_beliefs.dta", replace

*merge m:1 subject_id using "./Temp/ipclasses.dta" //risk aversion, blind prot choices and demographic vars
*drop _merge

*gen class_alt=2-class
*tab class class_alt
*label var class_alt "IP strategy class"
*label define clsnames 0 "Bayesians" 1 "Simpletons"
gen false_prob=phintWB+phintBW
label var false_prob "Prop. of lying gremlins"
*label values class_alt clsnames


set more off
eststo clear
eststo: reg bel_err phintBW phintWB i.subject_id, vce(cluster subject_id)
eststo: reg bel_err phintBW phintWB i.subject_id if blackhint==0, vce(cluster subject_id)
eststo: reg bel_err phintBW phintWB i.subject_id if blackhint==1, vce(cluster subject_id)
esttab using "./Tables/table_be_err.tex", b(%9.3f) se(%9.3f) ar2(%9.3f) label addnotes(Dep. variable: reported belief - posterior probability) title(Belief Elicitation: When Mistakes Happen) mtitles("All" "S=White" "S=Black") star("*" 0.10 "**" 0.05 "***" 0.01) indicate(Subject FE = *.subject_id) nobaselevels compress nogaps replace


*Testing FP/FN confusion:
eststo clear
eststo: reg be_ phintBW phintWB i.plevel i.subject_id, vce(cluster subject_id)
eststo: reg be_ phintBW phintWB i.plevel if blackhint==0, vce(cluster subject_id)
eststo: reg be_ phintBW phintWB i.plevel if blackhint==1, vce(cluster subject_id)
esttab using "./Tables/table_be_err1.tex", b(%9.3g) se(%9.1f) ar2(%9.2f) label addnotes(Dep. variable: reported belief - posterior probability) title(Belief Elicitation: When Mistakes Happen) mtitles("All" "S=White" "S=Black") star("*" 0.10 "**" 0.05 "***" 0.01) indicate(Subject FE = *.subject_id) nobaselevels compress nogaps replace


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
esttab using "./Tables/table_be3.tex", b(%9.3f) t(%9.1f) ar2(%9.2f) label drop(_cons) mtitles("OLS" "FE" "OLS" "FE" "OLS" "FE") star("*" 0.10 "**" 0.05 "***" 0.01) note("Decomposition works only for imperfect signals") nobaselevels compress nogaps replace

