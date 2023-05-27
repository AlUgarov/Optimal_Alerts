use "./Output/main_waves.dta", replace

keep subject_id participant_id round time_ip honest ncorrect p phintBW phintWB phintBB rev_response treatm_type tot_liars bl_gr w_gr
collapse (first) treatm_type participant_id time_ip honest ncorrect p phintBW phintWB phintBB rev_response tot_liars bl_gr w_gr, by(subject_id round)
label values treatm_type treatm_typel
save "./Temp/base_main_waves.dta", replace





use "./Output/main_waves.dta", replace

*Remove the timing information for now
drop *click* *page* history q104-q106 q114 q115
rename post_probW post_probw
rename post_probB post_probb

keep subject_id participant_id round ip_w ip_b be_w be_b be_time_w be_time_b post_probw post_probb time_ip honest ncorrect p phintBW phintWB phintBB rev_response


*Reshape to (subject_id round hint) long format
reshape long ip_ be_ be_time_ post_prob, i(subject_id round) j(hint) string

*Merging back subject-round variables:
merge m:1 subject_id round using "./Temp/base_main_waves.dta"
keep if _merge==3
drop _merge


**Merging risk aversion vars from the blind protection task
merge m:1 participant_id using "./Temp/blind_collapsed.dta"
keep if _merge==3
drop _merge

gen blackhint=.
replace blackhint=1 if hint=="b"
replace blackhint=0 if hint=="w"
gen question=2*(round-1)+blackhint

gen plevel=round(1000*p) //create integer var to code prior probability

sum(be_time_)
sum(time_ip)
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
gen accur_bel=tot_bel_err<=`err_med' //error is less than the median



*Evaluating relative belief accuracy conditional on tasks (prior x signal)
bys round: sum absbel_err
bys subject_id round: egen v1=sum(absbel_err)
bys p phintWB phintBW:egen med_bel_err=median(v1)
bys p phintWB phintBW:egen sd_bel_err=sd(v1)

gen accur_bel2=abs(bel_err)<=med_bel_err
gen accur_bel3=abs(bel_err)<0.005*sd_bel_err

sort  subject_id round hint
order subject_id round phintWB phintBW bel_err absbel_err accur_bel med_bel_err accur_bel2
bro


*Are subjects accurate in the second part of the experiment if they are accurate in the first?
gen laterounds=round>3
bys laterounds subject_id: egen sbel_err=sum(absbel_err)

bys laterounds: sum sbel_err




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
sum confid
bys plevel phintWB phintBW: sum confid


**Saving belief accuracy:
save "./Temp/base_main_waves.dta", replace


**ADDING HERE A BIT MORE ANALYSIS ON WHO IS MAKING THE MOST EGREGIOUS ERRORS**
use "./Temp/base_main_waves.dta", replace


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
sum minprotect
replace minprotect=1 if minprotect==2

gen post_proby=(1-ip_)*post_prob
by subject_id: egen maxrisk=max(post_proby)
sum maxrisk
gen nonmon=maxrisk>minprotect
tab nonmon
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


***There is little correlation btw making many mistakes on the quiz and the N of back switches in the BP!
tab goodquiz backswitcher0, chi2

*No correlation btw opposite errors in belief elicitation and back switches
tab backswitcher0 cert_err, chi2
tab backswitcher0 opp_err, chi2

*But there is a correlation btw quiz answers and large updating errors
tab goodquiz opp_err, chi2
reg cert_err goodquiz

*Obvious but imperfect correlation btw being a backswitcher and IP backswitcher!
tab backswitcher ip_backswitcher, chi2
tab backswitcher0 ip_backswitcher, chi2 //and here I'm doing the same for BP backswitchers without excluding repairable switches


*Obvious but imperfect correlation btw being a backswitcher and IP backswitcher!
tab goodquiz ip_backswitcher, chi2
tab goodquiz opp_ip_err, chi2 //and here I'm doing the same for BP backswitchers without excluding repairable switches


*weak correlation with backswitches:
tab backswitcher0 opp_ip_err, chi2
tab goodquiz opp_ip_err, chi2

gen badquiz=1-goodquiz
gen knuckleheadidness=badquiz+backswitcher0+ip_backswitcher+opp_err

use "./Temp/base_main_waves.dta", replace
collapse (mean) accur_bel accur_bel2 accur_bel3 be_change confid, by(subject_id round)
save "./Temp/bel_accuracy.dta", replace
