**BLIND PROTECTION ANALYSIS**
set more off
clear all


*!!put your working folder below:
cd C:\Tornado_warnings\Experiment\Alerts_Experiment
*cd C:\Tornado_warnings\Optimal_Alerts

set seed 135
//use "./Temp/allwaves_wide.dta", replace
use "./Output/demography.dta", replace

gen age23=age>23&!missing(age)
gen student=((educ==3)|(educ==4))&!missing(educ)
tab sex
tab stat_educ
tab student

hist correctcrt, width(`binwidth') fcolor(navy) lcolor(navy) title("Distribution of CRT scores") 
graph export "./Graphs/hist_crt.png", width(1200) height(800) replace


hist gpa, width(`binwidth') fcolor(navy) lcolor(navy) title("Distribution of GPA") 
graph export "./Graphs/hist_gpa.png", width(1200) height(800) replace


stop

file open mytbl using "Tables/demography.tex", write replace
file write mytbl "\begin{tabular}{l*{6}{c}}" _n
file write mytbl "\hline\hline" _n
file write mytbl " & \multicolumn{2}{|c|}{All}  & \multicolumn{2}{c|}{\$p \in\{0.1,0.3\}\$} & \multicolumn{2}{|c|}{\$p \in\{0.2,0.5\}\$}\\" _n "\hline" _n
file write mytbl " & N & \% & N & \% & N & \%  \\" _n "\hline" _n
file write mytbl "\multicolumn{7}{c}{All waves} \\" _n "\hline" _n
local rows sex age23 student stat_educ
recast double p
replace p=0.1 if abs(p-0.1)<0.01
replace p=0.2 if abs(p-0.2)<0.01
replace p=0.3 if abs(p-0.3)<0.01
replace p=0.5 if abs(p-0.5)<0.01
gen sample_13  = ((p==.1)|(p==.3))
gen sample_25  = ((p==.2)|(p==.5))
gen sample_all = (1==1)


//stop

label var sex "Male"
label var student "Students"
label var stat_educ "Had statistics classes"
label var age23 "Age\$>\$23yrs old"

foreach v of local rows {
    file write mytbl "`: var label `v'' "
    foreach s in all 13 25 {
        quietly count if `v' == 1 & sample_`s'
        local N = r(N)
        quietly count if sample_`s'
        local pct = 100 * `N' / r(N)
        file write mytbl "& `N' & `=round(`pct',1)' "
    }
    file write mytbl "\\" _n
}

file write mytbl "\hline" _n

file write mytbl "\multicolumn{7}{c}{First waves} \\" _n "\hline" _n


foreach v of local rows {
    file write mytbl "`: var label `v'' "
    foreach s in all 13 25 {
        quietly count if `v' == 1 & sample_`s' & wave==1
        local N = r(N)
        quietly count if sample_`s'
        local pct = 100 * `N' / r(N)
        file write mytbl "& `N' & `=round(`pct',1)' "
    }
    file write mytbl "\\" _n
}

file write mytbl "\hline" _n


file write mytbl "\multicolumn{7}{c}{Second wave} \\" _n "\hline" _n

foreach v of local rows {
    file write mytbl "`: var label `v'' "
    foreach s in all 13 25 {
        quietly count if `v' == 1 & sample_`s' & wave==2
        local N = r(N)
        quietly count if sample_`s'
        local pct = 100 * `N' / r(N)
        file write mytbl "& `N' & `=round(`pct',1)' "
    }
    file write mytbl "\\" _n
}



file write mytbl "\hline" _n "\end{tabular}"
file close mytbl
