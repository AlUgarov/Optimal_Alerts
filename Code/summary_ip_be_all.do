
*BASIC DATA ANALYSIS FOR BELIEF ELICITATION, AND INFORMED PROTECTION TASKS**
set more off
clear all

*!!put your working folder below:
cd C:\Tornado_warnings\Experiment\Alerts_Experiment
*cd C:\Tornado_warnings\Optimal_Alerts



*CHANGE TO THE PANEL STRUCTURE FOR THE OTHER TASKS********************
**Panel: participant_id round 
**because all the tasks except the blind protection have the same N of rounds (6) 
use "./Temp/secondwave_wide.dta", replace

*add the pilot's data:
*append using "./Temp/pilot_wide.dta"




reshape long ip_w_ ip_b_ ip_time_ be_w_ be_b_ be_time_w_ be_time_b_ wtp_1_ wtp_2_ wtp_3_ wtp_4_ wtp_5_ wtp_6_ wtp_7_ wtp_8_ wtp_9_ wtp_10_ wtp_11_ wtp_time_, i(participant_id sex age stat_educ gpa crt_payoff correctcrt) j(round)
rename *_ *
rename ip_time time_ip

**flipping priors orders for cases which were randomly flipped for new treatments in wave2:
gen chrono_round=round
recode round (1=4) (2=5) (3=6) (4=1) (5=2) (6=3) if (revorder==1)&!missing(revorder)
tab round
replace round=chrono_round





***Merging the treatment sequences***
merge m:1 seq round using "./Temp/exp_treatments_wave2.dta"
tab _merge
keep if _merge==3
drop _merge
drop if missing(participant_id)
gen wave=2
save "./Temp/secondwave_long.dta", replace

*stop


*Loading first waves data:
use "./Temp/mainwaves_wide.dta", replace

*add the pilot's data:
*append using "./Temp/pilot_wide.dta"
gen gpa=.
gen crt_payoff=.
gen correctcrt=.

reshape long ip_w_ ip_b_ ip_time_ be_w_ be_b_ be_time_w_ be_time_b_ wtp_1_ wtp_2_ wtp_3_ wtp_4_ wtp_5_ wtp_6_ wtp_7_ wtp_8_ wtp_9_ wtp_10_ wtp_11_ wtp_time_, i(participant_id sex age stat_educ gpa crt_payoff correctcrt) j(round)
rename *_ *
rename ip_time time_ip

***Merging the treatment sequences***
merge m:1 seq round using "./Temp/exp_treatments.dta"
drop _merge
gen wave=1

append using "./Temp/secondwave_long.dta"

encode participant_id, gen(subject_id) //as participant_id is initially a string

duplicates list participant_id round
duplicates list subject_id round




*stop

**Merging risk aversion vars from the blind protection task
merge m:1 participant_id using "./Temp/blind_collapsed_all.dta"
keep if _merge==3
drop _merge

drop if pilot==1 //dropping the pilot

sort participant_id round


save "./Output/all_waves.dta", replace

collapse (first) p sex age educ stat_educ fpayoff crt_payoff correctcrt gpa totprot switchprob wave, by(participant_id)

save "./Output/demography.dta", replace

use "./Output/all_waves.dta", replace




**Calculate the maximum willingness-to-pay for info for each participant/round (wtp):
foreach v of varlist wtp_1-wtp_11{
 display `v'
 replace `v'=0 if `v'==-99
}

egen wtp=rowtotal(wtp_1-wtp_11)
replace wtp=0.5*(wtp-1)
replace wtp=0 if wtp<0
sum wtp
*stop

gen honest_treatment=((bl_gr+w_gr)==0) //only honest gremlins in the group
replace be_w=0.01*(100-be_w) if pilot==0 //rescale belief elicitation responses to probabilities
replace be_w=0.01*be_w if pilot==1
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

replace ip_w=. if ip_w==-99
replace ip_b=. if ip_b==-99

gen ip_val=-(p*(phintWB*(1-ip_w)+phintBB*(1-ip_b))*loss+p*(phintWB*ip_w+phintBB*ip_b)*protectioncost+(1-p)*(phintWW*ip_w+phintBW*ip_b)*protectioncost)
replace ip_val=. if ip_b==-99
replace ip_val=. if ip_w==-99

sum ip_val



//Calculate the optimal protection strategy based on cost-loss ratio:
gen ip_w_o=0
replace ip_w_o=1 if post_probW>=protectioncost/loss
gen ip_b_o=0
replace ip_b_o=1 if post_probB>=protectioncost/loss


//Calculate exp. costs under the optimal strategy:
gen ip_val_o=-(p*(phintWB*(1-ip_w_o)+phintBB*(1-ip_b_o))*loss+p*(phintWB*ip_w_o+phintBB*ip_b_o)*protectioncost+(1-p)*(phintWW*ip_w_o+phintBW*ip_b_o)*protectioncost)

label var ip_val "Exp. costs"
label var ip_val_o "Optimal exp. costs"

gen ip_val_diff=ip_val_o-ip_val //discrepancy between actual and optimal expected costs of informed protection


//Calc exp costs under the optimal strategy for reported(!) beliefs:
gen phintB=p*phintBB+(1-p)*phintBW //ideally we should elicit these probabilities within the experiment
gen phintW=1-phintB

gen ip_w_mu=0
replace ip_w_mu=1 if be_w>=protectioncost/loss
gen ip_b_mu=0
replace ip_b_mu=1 if be_b>=protectioncost/loss

gen ip_val_mu=(loss*(phintW*be_w*(1-ip_w_mu)+phintB*be_b*(1-ip_b_mu))+protectioncost*(phintW*ip_w_mu+phintB*ip_b_mu))

gen rev_response=(ip_w==1)&(ip_b==0) //subjects protect when white and do not protect when black
tab rev_response


gen treatm_type=0
replace treatm_type=1 if bl_gr>0&w_gr==0
replace treatm_type=2 if w_gr>0&bl_gr==0
replace treatm_type=3 if bl_gr>0&w_gr>0
label define treatm_typel 0 "Honest" 1 "White-eyed only" 2 "Black-eyed only" 3 "All types"
label values treatm_type treatm_typel

*clean out the subjects with missing informed protection responses (1)
gen miss_ip=missing(ip_w)|missing(ip_b)
bys subject_id: egen incompl_ip=sum(miss_ip)
drop if incompl_ip>0



duplicates list subject_id round //no duplicates

tab wave

*Saving the cleaned dataset with the panel structure
save "./Output/all_waves.dta", replace







use "./Output/all_waves.dta", replace
*stop
**Create the summary statistics table:
gen duration_min=duration/60
collapse (first) seq sex age stat_educ college final_payoff duration_min, by(subject_id)

gen seq_type=(seq>3)

gen old=age>23
replace college=1-college

tabstat sex old college stat_educ, by(seq_type) statistics(sum mean) column(statistics) longstub


use "./Output/second_wave_test.dta", replace

keep subject_id participant_id round time_ip honest ncorrect p phintBW phintWB phintBB rev_response treatm_type tot_liars bl_gr w_gr
collapse (first) treatm_type participant_id time_ip honest ncorrect p phintBW phintWB phintBB rev_response tot_liars bl_gr w_gr, by(subject_id round)
label values treatm_type treatm_typel

duplicates list subject_id round //no duplicates


save "./Temp/base_second_wave.dta", replace



use "./Output/second_wave_test.dta", replace
*Remove the timing information for now
drop *click* *page* history q104-q106 q114 q115
rename post_probW post_probw
rename post_probB post_probb

keep subject_id participant_id round ip_w ip_b be_w be_b be_time_w be_time_b post_probw post_probb time_ip honest ncorrect p phintBW phintWB phintBB rev_response wave


duplicates list participant_id round //no duplicates


*replace subject_id=subject_id+100000*(wave-1)

*Reshape to (subject_id round hint) long format
reshape long ip_ be_ be_time_ post_prob, i(subject_id round) j(hint) string
*stop

*Merging back subject-round variables:
merge m:1 subject_id round using "./Temp/base_second_wave.dta"
keep if _merge==3
drop _merge


duplicates list participant_id round hint
*stop


**Merging risk aversion vars from the blind protection task
merge m:1 participant_id using "./Temp/blind_collapsed_all.dta"
keep if _merge==3
drop _merge

gen blackhint=.
replace blackhint=1 if hint=="b"
replace blackhint=0 if hint=="w"
gen question=2*(round-1)+blackhint

gen plevel=round(1000*p) //create integer var to code prior probability

label var ip_ "Informed protection"
label var be_ "Belief"
label var blackhint "S=Black"
label var honest "Honest treatment"
label var post_prob "Posterior prob."

replace ip_=. if ip_==-99


**Calculate subject-level accuracy of reported beliefs:
gen bel_err=post_prob-be_
label var bel_err "Belief error"

gen absbel_err=abs(bel_err)
sort subject_id
by subject_id: egen tot_bel_err=sum(absbel_err) //total abs error per subject

sum tot_bel_err, detail
local err_med=r(p50)
gen accur_bel=tot_bel_err<`err_med' //error is less than the median



*Evaluating relative belief accuracy conditional on tasks (prior x signal)
bys round: sum absbel_err
bys subject_id round: egen v1=sum(absbel_err)
bys p phintWB phintBW:egen med_bel_err=median(v1)
bys p phintWB phintBW:egen sd_bel_err=sd(v1)
gen accur_bel2=abs(bel_err)<med_bel_err
gen accur_bel3=abs(bel_err)<0.1*sd_bel_err

sort  subject_id round hint
order subject_id round phintWB phintBW bel_err absbel_err accur_bel med_bel_err accur_bel2


*Are subjects accurate in the second part of the experiment if they are accurate in the first?
gen laterounds=round>3
bys laterounds subject_id: egen sbel_err=sum(absbel_err)

*How much the signal affects beliefs?
gen be_w=be_
gen be_b=be_
replace be_w=. if hint=="b"
replace be_b=. if hint=="w"
bys subject_id round: egen be_w1=mean(be_w)
bys subject_id round: egen be_b1=mean(be_b)
replace be_w=be_w1
replace be_b=be_b1

gen be_change=be_b-be_w



*Calculate average confidence as reported beliefs that the decision is correct for each signal:
gen post_prob_w=1-post_prob
gen post_prob_b=post_prob
replace post_prob_w=. if hint=="b"
replace post_prob_b=. if hint=="w"
bys subject_id round: egen post_prob_w1=mean(post_prob_w)
bys subject_id round: egen post_prob_b1=mean(post_prob_b)

gen prob_correct=be_
replace prob_correct=1-be_ if hint=="w"
gen prob_B_sign=(p*phintBB+(1-p)*phintBW)
gen prob_W_sign=1-prob_B_sign

gen confid=prob_W_sign*(1-be_w1)+prob_B_sign*be_b1

duplicates list subject_id round hint
*stop

**Saving belief accuracy:
save "./Temp/base_second_wave.dta", replace
*stop

**ADDING HERE A BIT MORE ANALYSIS ON WHO IS MAKING THE MOST EGREGIOUS ERRORS**
use "./Temp/base_second_wave.dta", replace


drop if abs(0.5-post_prob)<0.499 //keeping only certain signals

gen opp_err0=abs(be_-post_prob)>0.999 //giving the opposite probability (0 vs 1)
gen cert_err0=(be_>0)&(abs(0.5-post_prob)>0.499) //belief error when the outcome is completely certain cond. on a signal

bys subject_id: egen cert_err=sum(cert_err0)
bys subject_id: egen opp_err=max(opp_err0)



**IP consistency: non-monotonicies, protecting when p=0, not protecting when p=1
sort subject_id post_prob

*calculate the posterior at the earliest protection
gen post_probx=post_prob
replace post_probx=2 if ip_==1
by subject_id: egen minprotect=min(post_probx)
replace minprotect=1 if minprotect==2

gen post_proby=(1-ip_)*post_prob
by subject_id: egen maxrisk=max(post_proby)
sum maxrisk
gen nonmon=maxrisk>minprotect

*then find if there are no protection for any higher probabilities
gen backswitch=(ip_==0&ip[_n-1]==1)&(subject_id==subject_id[_n-1])

*How to check more properly if it's monotonic?
tab backswitch
bys subject_id: egen nip_backswitches=sum(backswitch)
tab nip_backswitches
bys subject_id: egen ip_backswitcher=max(backswitch)


gen toosafe=(ip_==1)&(post_prob==0)
gen suicidal=(ip_==0)&(post_prob==1)
gen opp_ip_err0=toosafe+suicidal
by subject_id: egen opp_ip_err=max(opp_ip_err0)
tab opp_ip_err

bys laterounds: egen opp_ip_err1=max(opp_ip_err0)
tab opp_ip_err1

*Comparing early rounds with later rounds (0.3 and 0.5 priors):
bys laterounds: sum opp_ip_err0 //IP

bys laterounds: sum cert_err0 //BE when the posterior should be certain
bys laterounds: sum opp_err0 //Reporting the opposite belief when the posterior is certain, slightly lower


collapse (first) backswitcher backswitcher0 nbswitches goodquiz opp_ip_err ip_backswitcher opp_err cert_err, by(subject_id)

gen badquiz=1-goodquiz
gen knuckleheadidness=badquiz+backswitcher0+ip_backswitcher+opp_err

*The histogram doesn't indicate any cluster of "knuckleheads", but most subjects make at least one large mistake
hist knuckleheadidness, start(0) discrete frequency title("Subjects making N large errors") note("Counted mistakes: 2 or more quiz mistakes, back switches in IP or BP," " stating opposite beliefs when the color is certain.") color(navy)
graph export "./Graphs/hist_bigerrors_test.png", width(1200) height(800) replace



use "./Temp/base_second_wave.dta", replace
collapse (mean) accur_bel accur_bel2 accur_bel3 be_change confid absbel_err tot_bel_err, by(subject_id round)

save "./Temp/bel_accuracy_all.dta", replace

use "./Temp/base_second_wave.dta", replace


label define accur_bel_l 0 "Inaccurate" 1 "Accur. beliefs"
label values accur_bel accur_bel_l
label values accur_bel2 accur_bel_l

label values stat_educ stat_educl

gen small_error=abs(bel_err)<0.1

gen certaincolor=abs(0.5-post_prob)>0.499

sum absbel_err if abs(0.5-post_prob)<0.499, detail

**Exploring belief elicitation errors for both low and high priors:
hist bel_err if p<0.299, title("Errors in elicited beliefs (low priors)") xtitle("Posterior - Belief") fraction note("By belief elicitation task, no aggregation to round or subjects") color(navy)
graph export "./Graphs/hist_belief_error_low_w2.png", width(1200) height(800) replace

hist bel_err if p>0.299, title("Errors in elicited beliefs (high priors)") xtitle("Posterior - Belief") fraction note("By belief elicitation task, no aggregation to round or subjects") color(navy)
graph export "./Graphs/hist_belief_error_high_w2.png", width(1200) height(800) replace



**Do a more formal testing of equality of variances btw low and high priors
bys laterounds: sum bel_err
sdtest bel_err, by(laterounds)
robvar bel_err, by(laterounds)

**Is it just because priors are higher? Test it for groups with low priors vs high priors:
gen seqtype=(plevel==200)|(plevel==500)
tab seqtype
sdtest bel_err, by(seqtype)  //the difference in variances is not significant
*stop


hist bel_err, title("Errors in elicited beliefs") xtitle("Posterior - Belief") fraction note("By belief elicitation task, no aggregation to round or subjects") color(navy)
graph export "./Graphs/hist_belief_error.png", width(1200) height(800) replace

hist bel_err if abs(0.5-post_prob)<0.499, title("Errors in beliefs, ball color is uncertain") xtitle("Posterior - Belief") note("By belief elicitation task, no aggregation to round or subjects") fraction color(navy)
graph export "./Graphs/hist_belief_error_s4.png", width(1200) height(800) replace

hist bel_err if abs(0.5-post_prob)>0.499, title("Errors in beliefs, ball color is certain") xtitle("Posterior - Belief") note("By belief elicitation task, no aggregation to round or subjects") fraction color(navy)
graph export "./Graphs/hist_belief_error_s5.png", width(1200) height(800) replace



qui reg be_ post_prob
local r2 : display %5.3f = e(r2)
corr be_ post_prob
local rho:  display %5.3f = r(rho)
graph twoway (scatter be_ post_prob, jitter(1)) (lfit be_ post_prob) , title("Belief updating") xtitle("True probability") ytitle("Elicited belief")  note("All obs, correlation=`rho'") legend(off)
graph export "./Graphs/updating_s1.png", width(1200) height(800) replace

qui reg be_ post_prob if goodquiz==1
local r2 : display %5.3f = e(r2)
corr be_ post_prob if pilot==0&goodquiz==1
local rho:  display %5.3f = r(rho)
graph twoway  (scatter be_ post_prob, jitter(1)) (lfit be_ post_prob) if pilot==0&goodquiz==1, title("Belief updating") xtitle("True probability") ytitle("Elicited belief") legend(off) note("Good quiz, correlation=`rho'")
graph export "./Graphs/updating_s2.png", width(1200) height(800) replace

qui reg be_ post_prob if honest==0
local r2 : display %5.3f = e(r2)
corr be_ post_prob if pilot==0&honest==0
local rho:  display %5.3f = r(rho)
graph twoway  (scatter be_ post_prob, jitter(1)) (lfit be_ post_prob) if pilot==0&honest==0, title("Belief updating") xtitle("True probability") ytitle("Elicited belief") legend(off) note("Excluding certain signals, correlation=`rho'")
graph export "./Graphs/updating_s3.png", width(1200) height(800) replace

qui reg be_ post_prob if abs(0.5-post_prob)<0.499
local r2 : display %5.3f = e(r2)
corr be_ post_prob if pilot==0&abs(0.5-post_prob)<0.499
local rho:  display %5.3f = r(rho)
graph twoway (scatter be_ post_prob, jitter(1)) (lfit be_ post_prob) if pilot==0&abs(0.5-post_prob)<0.499, title("Belief vs posterior, ball color is uncertain") xtitle("True probability") ytitle("Elicited belief") legend(off) note("Correlation=`rho'")
graph export "./Graphs/updating_s4.png", width(1200) height(800) replace

qui reg be_ post_prob if abs(0.5-post_prob)>0.499
local r2 : display %5.3f = e(r2)
corr be_ post_prob if pilot==0&abs(0.5-post_prob)>0.499
local rho:  display %5.3f = r(rho)
graph twoway  (scatter be_ post_prob, jitter(1)) (lfit be_ post_prob) if pilot==0&abs(0.5-post_prob)>0.499, title("Belief vs posterior, ball color is certain") xtitle("True probability") ytitle("Elicited belief") legend(off) note("Correlation=`rho'")
graph export "./Graphs/updating_s5.png", width(1200) height(800) replace

save "./Temp/base_main_waves.dta", replace


