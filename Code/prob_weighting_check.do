

*BASIC DATA ANALYSIS FOR BELIEF ELICITATION, AND INFORMED PROTECTION TASKS**
set more off
clear all

capture program drop gw_dphi

program define gw_dphi

    version 18.0

    /*
    Creates phi'(x) for the Gonzalez-Wu weighting function:

                         delta*gamma*x^(gamma-1)*(1-x)^(gamma-1)
        phi'(x) = ---------------------------------------------------
                   [delta*x^gamma + (1-x)^gamma]^2
    */

    syntax newvarname, ///
        AT(varname numeric) ///
        DELta(varname numeric) ///
        GAMma(varname numeric)

    quietly {

        assert `at' > 0 & `at' < 1 if !missing(`at')
        assert `delta' > 0 if !missing(`delta')
        assert `gamma' > 0 if !missing(`gamma')

        generate double `varlist' =                         ///
            `delta' * `gamma'                              ///
            * `at'^(`gamma' - 1)                           ///
            * (1 - `at')^(`gamma' - 1)                    ///
            / (`delta' * `at'^`gamma'                      ///
               + (1 - `at')^`gamma')^2
    }

end

*!!put your working folder below:*
*cd C:\Tornado_warnings\Experiment\Alerts_Experiment
*cd C:\Tornado_warnings\Optimal_Alerts
cd "/Users/alexander/Optimal_Alerts"

import delimited using "./Input/gonzalez_wu_estimates.csv", varnames(1) clear

gen id=1
*will use much later to loop over:
levelsof subject_id, local(gw_subjects)
reshape wide alpha delta gamma, i(id) j(subject_id)

save "./Temp/gonzalez_wu_estimates.dta", replace



use "./Temp/all_treatments.dta", replace

*the file has phintWB as probability of false-negative (P_{01}) and phintBW as the false-positive probability (P_{10}), p is \pi*

gen id=1
merge m:1 id using "./Temp/gonzalez_wu_estimates.dta"

assert _merge==3
drop _merge

gen phintWW = 1 - phintBW

* Probability of the best outcome:
* Pr(Y-b) = (1-pi)P00
gen cum_best = (1-p)*phintWW

* Cumulative probability of the best and middle outcomes:
* Pr(Y-b or Y-b-c) = 1-pi*P01
gen cum_best_middle = 1 - p*phintWB

label variable cum_best ///
    "(1-p)P00: cumulative probability of best outcome"

label variable cum_best_middle ///
    "1-p*P01: cumulative probability of best and middle outcomes"
	
	
gen byte pw_case=0



* Ordinary interior calculation
replace pw_case=1 if cum_best > 0 & cum_best_middle < 1
gen n_lambda_gt1=0
foreach s of local gw_subjects {
    di `s'
    gw_dphi dphi_fn`s' if pw_case == 1, at(cum_best_middle) delta(delta`s') gamma(gamma`s')

    gw_dphi dphi_fp`s' if pw_case == 1, at(cum_best) delta(delta`s') gamma(gamma`s')

    gen double lambda`s' = dphi_fn`s'/dphi_fp`s' if pw_case == 1

    gen byte lambda_gt1_`s' = lambda`s' > 1 if pw_case == 1
	replace n_lambda_gt1=n_lambda_gt1+lambda_gt_`s'
}

	
sum n_lambda_gt1 



