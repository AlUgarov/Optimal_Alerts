
**BLIND PROTECTION ANALYSIS**
set more off
clear all


*!!put your working folder below:
*cd C:\Tornado_warnings\Experiment\Alerts_Experiment
*cd C:\Tornado_warnings\Optimal_Alerts

set seed 135

use "./Temp/mainwaves_wide.dta", replace
* bp - protection decision (0 - do not protect, 1 - protect)

keep participant_id bp_* bp_time_* sex age stat_educ ncorrect informed_correct wtpq_correct educ
reshape long bp_ bp_time_,  i(participant_id sex age stat_educ ncorrect) j(round)
rename bp_ bp
rename bp_time_ submittime
gen p=0.05*round


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

*if only one back switch and it can be repaired by a single change
*then change the switchround to the 7-total number of round using protection
gen repairable=0
replace repairable=1 if (nbswitches==1)&(bp[_n-1]!=bp)&(bp[_n+1]!=bp)&(round>1)&(round<6)
replace repairable=1 if (nbswitches==1)&(bp[_n-1]!=bp)&(bp==0)&(round==6)


save "./Temp/bp_val.dta", replace



gen college=educ>4

*Collapsing to have one obs per participant (study participant's characteristics):
collapse (mean) bp submittime (first) sex age stat_educ ncorrect informed_correct wtpq_correct college (sum) totprot=bp (max) switcher backswitcher nbswitches repairable (min) maxspeed=submittime firstswitch=switchround backswitchround, by(participant_id)

tab firstswitch
replace firstswitch=. if backswitcher==1
replace firstswitch=7-totprot if repairable==1
tab firstswitch

gen backswitcher0=backswitcher
replace backswitcher=0 if repairable==1
gen byte allprotect=(totprot==6)


**Estimate subjects' risk aversion from the blind protection choices:
gen switchprob=0.05*firstswitch
replace switchprob=-0.5 if missing(switchprob)
gen pilot=0



**Merge-in blind protection choices from the pilot
*append using "./Temp/blind_collapsed_pilot.dta"

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


save "./Temp/blind_collapsed.dta", replace

use "./Temp/bp_val.dta", replace
gen pilot=0

*append using "./Temp/bp_val_pilot.dta"
replace pilot=1 if missing(pilot)

keep participant_id p bp bp_val pilot

save "./Temp/bp_val.dta", replace
gen xid=1
collapse (mean) mean=bp (count) n=xid, by(p)
gen se=sqrt(mean*(1-mean)/n)


*Blind protection average response diagram (Fig. 1)
serrbar mean se p, scale (1.96) title("Blind Protection Response") lwidth(thick) xtitle("Probability of a black ball") ytitle("Proportion of protection choices") ysc(r(0 1)) ylabel(0(0.2)1.0)
graph export "./Graphs/blind_prot_sta.png", width(1000) height(1000) replace

gen task_type="Blind"
rename p post_prob
keep mean post_prob se task_type
save "./Temp/bp_graph_data.dta", replace
