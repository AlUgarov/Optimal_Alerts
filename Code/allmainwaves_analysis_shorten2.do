
**ANALYZE THE FIRST WAVE (planned 100 participants)***
set more off
clear all

**Requires: estout, moremata

*!!put your working folder below:
cd C:\Tornado_warnings\Experiment\Alerts_Experiment
*cd C:\Tornado_warnings\Optimal_Alerts

set seed 135

*Prepare the treatment characteristics to merge:
import delimited using "./Input/exp_treatments_pilot.csv", varnames(1) clear //I prepare this file separately in R, describes each potential treatment (prior prob+signal information structure)
drop v1
rename snames treatn

save "./Temp/exp_treatments_pilot.dta", replace

import delimited using "./Input/treatment_sequences.csv", varnames(1) clear //I prepare this file separately in Excel, describes treatment order for each group of subjects/session type
rename ïseq seq
reshape long r, i(seq) j(round) //so it is sequence - round structure to make merging with the results file easier
rename r treatn

merge m:1 treatn using "./Temp/exp_treatments_pilot.dta" //merging characteristics of each treatment (prior prob, info structure = N of each gremlin types)
drop _merge
sort seq round


*Create additional variables to describe posterior probabilities:
gen phintWB=(w_gr/tot_gr)
gen phintWW=(tot_gr-bl_gr)/tot_gr
gen phintBB=1-phintWB
gen phintBW=1-phintWW
gen post_probW=p*phintWB/(p*phintWB+(1-p)*phintWW) //posterior prob that the ball is black if the HINT is white
gen post_probB=p*phintBB/(p*phintBB+(1-p)*phintBW) //posterior prob that the ball is black if the HINT is black

gen tot_liars=bl_gr+w_gr

save "./Temp/exp_treatments.dta", replace


*First importing the Qualtrics data
*The results data has variable names incompatible with Stata
*Hence we import the first row of the results file first as data. The first row contains original variable names, which we process to an acceptable format and store in a string 
import delimited using "./Input/Data_All_mainwaves.csv", rowrange(1:1) varnames(nonames) clear
local myvarlist
foreach v of varlist _all{
 display `v'
 replace `v'=strtrim(`v')
 replace `v'=subinstr(`v', " (in seconds)", "", .)
 replace `v'=subinstr(`v', " ", "", .)
 replace `v'=substr(`v',3,.)+"_"+substr(`v',1,1) if regexm(`v',"^[0-9]")
 replace `v'=strlower(`v')
 local luck=`v'[1]
 local myvarlist `myvarlist' `luck'
}
*Then we import the data from the results file (ignoring first 3 rows with var names and definition). Then we rename these variables using the previously stored processed names 
import delimited using "./Input/Data_All_mainwaves.csv", rowrange(4:110) varnames(nonames) clear

display "`myvarlist'"
local i=1
foreach v of varlist _all{
 di `i'
 local luck: word `i' of `myvarlist'
 di "`luck'"
 rename `v' `luck'
 local i=`i'+1
}


*Mass renaming some response variables to make reshaping and analysis more convenient:
rename b1_* bp_* //blind protection responses by round
rename be1_b_1* be_b* //beliefs if the gremlin says "black" by round
rename be1_w_15* be_w* //beliefs if the gremlin says "white" by round
rename wtp1_*_1_* wtp_*_* //value responses by round

rename q59_pagesubmit_* bp_time_* //blind protection time to submit
rename q61_pagesubmit_* ip_time_* //informed protection time to submit
rename q79_pagesubmit_* be_time_w_* //belief elicitation time to submit (white hint)
rename q64_pagesubmit_* be_time_b_* //belief elicitation time to submit (black hint)
rename q185_pagesubmit_* wtp_time_* //wtp elicitation time to submit
rename sequence seq




***CLEANING************************
drop if missing(participant_id) //if no id
keep if finished==1 //subjects complete the experiment
drop ipaddress recipient* distributionchannel userlanguage finished //drop unnecessary identifying information
drop round externalreference location* responseid t1r1 sel* treatment_plans ip_wc ip_bc payoffs signal* //other excessive info, relicts of previous experiment edits

recode q6 (1=0) (2=1) (4=0), gen(sex)
label var sex "Male"
label define sexes 0 "Non-male" 1 "Male"
label values sex sexes

gen age=2021-(2003-75)-q8 //18 means 18 and younger!!

recode q137 (1=1) (2=0) (3=0),gen(stat_educ)
label var stat_educ "Stat. class"
label define stat_educl 0 "No" 1 "Stat. class"
label values stat_educ stat_educl

rename q10 educ
label var educ "Education level"
label define educl 1 "No schooling" 2 "Grades 1-12, no diploma" 3 "High school diploma or GED" 4 "Some college, but no degree" 5 "Associate or bachelor's degree" 6 "Graduate or professional degree"
label values educ educl


*identify correct quiz answers
gen blind_correct=(q157==2)+(q158==3)+(q189==4)
gen informed_correct=(q120==2)+(q119==2)+(q121==1)+(q118==3)+(q117==4)
gen add_inf_corect=(q130==1)+(q135==3)+(q136==3)
gen wtpq_correct=(q164==3)+(q166==3)
gen ncorrect=blind_correct+informed_correct+wtpq_correct

gen pilot=0

save "./Temp/mainwaves_wide.dta", replace



import delimited using "./Input/coded_strategies.csv", clear

gen mention_hint=0
replace mention_hint=1 if strpos(explanation, "hint")>0
replace mention_hint=1 if strpos(explanation, "says")>0
replace mention_hint=1 if strpos(explanation, "say")>0
replace mention_hint=1 if strpos(explanation, "tell")>0
replace mention_hint=1 if strpos(explanation, "told")>0

gen mention_prop=0
replace mention_prop=1 if strpos(explanation, "balls")

gen mention_hon=0
replace mention_hon=1 if strpos(explanation, "truth")>0
replace mention_hon=1 if strpos(explanation, "honest")>0
replace mention_hon=1 if strpos(explanation, "lying")>0


save "./Temp/coded_strategies.dta", replace
use "./Temp/mainwaves_wide.dta", replace
merge 1:1 participant_id using "./Temp/coded_strategies.dta"
drop _merge



save "./Temp/mainwaves_wide.dta", replace

**BLIND PROTECTION ANALYSIS**
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
append using "./Temp/blind_collapsed_pilot.dta"

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
serrbar mean se p, scale (1.96) title("Blind Protection Response") lwidth(thick) xtitle("Probability of a black ball") ytitle("Proportion of protection choices") ysc(r(0 1)) ylabel(0(0.2)1.0) note("The bars show 95% confidence intervals for the mean proportion of subjects " "choosing protection at each probability.")
graph export "./Graphs/blind_prot_sta.png", width(1000) height(1000) replace

gen task_type="Blind"
rename p post_prob
keep mean post_prob se task_type
save "./Temp/bp_graph_data.dta", replace




**CHANGE TO THE PANEL STRUCTURE FOR THE OTHER TASKS********************
**Panel: participant_id round 
**because all the tasks except the blind protection have the same N of rounds (6) 
use "./Temp/mainwaves_wide.dta", replace

*add the pilot's data:
*append using "./Temp/pilot_wide.dta"
encode participant_id, gen(subject_id) //as participant_id is initially a string

cluster kmeans ip_w_1 ip_b_1 ip_w_2 ip_b_2 ip_w_3 ip_b_3, k(3)



reshape long ip_w_ ip_b_ ip_time_ be_w_ be_b_ be_time_w_ be_time_b_ wtp_1_ wtp_2_ wtp_3_ wtp_4_ wtp_5_ wtp_6_ wtp_7_ wtp_8_ wtp_9_ wtp_10_ wtp_11_ wtp_time_, i(subject_id participant_id sex age stat_educ) j(round)
rename *_ *
rename ip_time time_ip

***Merging the treatment sequences***
merge m:1 seq round using "./Temp/exp_treatments.dta"
drop _merge


**Merging risk aversion vars from the blind protection task
merge m:1 participant_id using "./Temp/blind_collapsed.dta"
keep if _merge==3
drop _merge

drop if pilot==1 //dropping the pilot

sort participant_id round

**Calculate the maximum willingness-to-pay for info for each participant/round (wtp):
foreach v of varlist wtp_1-wtp_11{
 display `v'
 replace `v'=0 if `v'==-99
}

egen wtp=rowtotal(wtp_1-wtp_11)
replace wtp=0.5*(wtp-1)
replace wtp=0 if wtp<0
sum wtp

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

gen ip_val_diff=ip_val-ip_val_o //discrepancy between actual and optimal expected costs of informed protection

replace ip_val_diff=-ip_val_diff


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

*Saving the cleaned dataset with the panel structure
save "./Output/main_waves.dta", replace
use "./Output/main_waves.dta", replace

**Create the summary statistics table:
gen duration_min=duration/60
collapse (first) seq sex age stat_educ college final_payoff duration_min, by(subject_id)

gen seq_type=(seq>3)

gen old=age>23
replace college=1-college

tabstat sex old college stat_educ, by(seq_type) statistics(sum mean) column(statistics) longstub
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
hist tot_bel_err
gen accur_bel=tot_bel_err<`err_med' //error is less than the median



*Evaluating relative belief accuracy conditional on tasks (prior x signal)
bys round: sum absbel_err
bys subject_id round: egen v1=sum(absbel_err)
bys p phintWB phintBW:egen med_bel_err=median(v1)
gen accur_bel2=abs(bel_err)<med_bel_err

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

*The histogram doesn't indicate any cluster of "knuckleheads", but most subjects make at least one large mistake
hist knuckleheadidness, start(0) discrete frequency title("Subjects making N large errors") note("Counted mistakes: 2 or more quiz mistakes, back switches in IP or BP," " stating opposite beliefs when the color is certain.") color(navy)
graph export "./Graphs/hist_bigerrors.png", width(1200) height(800) replace



use "./Temp/base_main_waves.dta", replace
collapse (mean) accur_bel accur_bel2 be_change confid, by(subject_id round)
save "./Temp/bel_accuracy.dta", replace
use "./Temp/base_main_waves.dta", replace



label define accur_bel_l 0 "Inaccurate" 1 "Accur. beliefs"
label values accur_bel accur_bel_l
label values accur_bel2 accur_bel_l

label values stat_educ stat_educl

gen small_error=abs(bel_err)<0.1

gen certaincolor=abs(0.5-post_prob)>0.499

sum absbel_err if abs(0.5-post_prob)<0.499, detail

**Exploring belief elicitation errors for both low and high priors:
hist bel_err if p<0.299, title("Errors in elicited beliefs (low priors)") xtitle("Posterior - Belief") fraction note("By belief elicitation task, no aggregation to round or subjects") color(navy)
graph export "./Graphs/hist_belief_error_low.png", width(1200) height(800) replace

hist bel_err if p>0.299, title("Errors in elicited beliefs (high priors)") xtitle("Posterior - Belief") fraction note("By belief elicitation task, no aggregation to round or subjects") color(navy)
graph export "./Graphs/hist_belief_error_high.png", width(1200) height(800) replace



**Do a more formal testing of equality of variances btw low and high priors
bys laterounds: sum bel_err
sdtest bel_err, by(laterounds)
robvar bel_err, by(laterounds)

**Is it just because priors are higher? Test it for groups with low priors vs high priors:
gen seqtype=(plevel==200)|(plevel==500)
tab seqtype
sdtest bel_err, by(seqtype)  //the difference in variances is not significant


hist bel_err, title("Errors in elicited beliefs") xtitle("Posterior - Belief") fraction note("By belief elicitation task, no aggregation to round or subjects") color(navy)
graph export "./Graphs/hist_belief_error.png", width(1200) height(800) replace

hist bel_err  if abs(0.5-post_prob)<0.499, title("Errors in elicited beliefs") xtitle("Posterior - Belief") fraction note("Main waves only, excluding certain signals") color(navy)
graph export "./Graphs/hist_belief_error_s3.png", width(1200) height(800) replace

hist bel_err  if abs(0.5-post_prob)<0.499, title("Errors in beliefs, ball color is uncertain") xtitle("Posterior - Belief") fraction note("Main waves only") color(navy)
graph export "./Graphs/hist_belief_error_s4.png", width(1200) height(800) replace

hist bel_err  if abs(0.5-post_prob)>0.499, title("Errors in beliefs, ball color is certain") xtitle("Posterior - Belief") fraction note("Main waves only") color(navy)
graph export "./Graphs/hist_belief_error_s5.png", width(1200) height(800) replace

qui reg be_ post_prob
local r2 : display %5.3f = e(r2)
corr be_ post_prob
local rho:  display %5.3f = r(rho)
graph twoway (scatter be_ post_prob, jitter(1)) (lfit be_ post_prob) , title("Belief updating") xtitle("True probability") ytitle("Elicited belief")  note("All obs including pilot, correlation=`rho'") legend(off)
graph export "./Graphs/updating_s1.png", width(1200) height(800) replace

qui reg be_ post_prob if pilot==0&goodquiz==1
local r2 : display %5.3f = e(r2)
corr be_ post_prob if pilot==0&goodquiz==1
local rho:  display %5.3f = r(rho)
graph twoway  (scatter be_ post_prob, jitter(1)) (lfit be_ post_prob) if pilot==0&goodquiz==1, title("Belief updating") xtitle("True probability") ytitle("Elicited belief") legend(off) note("Main waves only, good quiz, correlation=`rho'")
graph export "./Graphs/updating_s2.png", width(1200) height(800) replace

qui reg be_ post_prob if pilot==0&honest==0
local r2 : display %5.3f = e(r2)
corr be_ post_prob if pilot==0&honest==0
local rho:  display %5.3f = r(rho)
graph twoway  (scatter be_ post_prob, jitter(1)) (lfit be_ post_prob) if pilot==0&honest==0, title("Belief updating") xtitle("True probability") ytitle("Elicited belief") legend(off) note("Main waves only, excluding certain signals, correlation=`rho'")
graph export "./Graphs/updating_s3.png", width(1200) height(800) replace

qui reg be_ post_prob if pilot==0&abs(0.5-post_prob)<0.499
local r2 : display %5.3f = e(r2)
corr be_ post_prob if pilot==0&abs(0.5-post_prob)<0.499
local rho:  display %5.3f = r(rho)
graph twoway (scatter be_ post_prob, jitter(1)) (lfit be_ post_prob) if pilot==0&abs(0.5-post_prob)<0.499, title("Belief vs posterior, ball color is uncertain") xtitle("True probability") ytitle("Elicited belief") legend(off) note("Main waves only, correlation=`rho'")
graph export "./Graphs/updating_s4.png", width(1200) height(800) replace

qui reg be_ post_prob if pilot==0&abs(0.5-post_prob)>0.499
local r2 : display %5.3f = e(r2)
corr be_ post_prob if pilot==0&abs(0.5-post_prob)>0.499
local rho:  display %5.3f = r(rho)
graph twoway  (scatter be_ post_prob, jitter(1)) (lfit be_ post_prob) if pilot==0&abs(0.5-post_prob)>0.499, title("Belief vs posterior, ball color is certain") xtitle("True probability") ytitle("Elicited belief") legend(off) note("Main waves only, correlation=`rho'")
graph export "./Graphs/updating_s5.png", width(1200) height(800) replace


**********************************************
****----MAIN REGRESSIONS (beliefs and informed protection)---***********
**********************************************
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
sum tot_diff
list subject_id totdisagr1 totdisagr2 ip_dev tot_diff if ip_==0&post_prob==1

*Notes: nonsensical responses mostly come from subjects choosing different answers to different probabilities, 
* most subjects choosing no protection when prob=1 do not choose protection when prob=0
* seems to be that these subjects are just more random: sum of deviations btw posteriors and responses is much higher for them than for others even excluding that one mistake

collapse (mean) mean=ip_ (count) n=xid, by(post_probi)
gen se=sqrt(mean*(1-mean)/n)

gen task_type="Informed"



gen post_prob=0.001*post_probi
append using "./Temp/bp_graph_data.dta"
serrbar mean se post_prob if task_type=="Informed", scale (1.96) title("Informed Protection Response") lwidth(thick) xtitle("Posterior probability of a black ball") ytitle("Proportion of protection choices")  note("The bars show 95% confidence intervals for the mean proportion of subjects " "choosing protection at each probability.")
graph export "./Graphs/ip_response.png", width(1000) height(1000) replace



gen low=mean-1.96*se
gen high=mean+1.96*se
twoway (rcap high low post_prob if task_type=="Blind", lwidth(thick)) (rcap high low post_prob if task_type=="Informed"), title("Protection Response") xtitle("Posterior probability of a black ball") ytitle("Proportion of protection choices") legend(label(1 "Blind") label(2 "Informed"))
graph export "./Graphs/ip_response_comp.png", width(1200) height(800) replace

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
gen fef=1
eststo clear
eststo: logit ip_ phintBW phintWB highprob, vce(cluster subject_id)
local r2p=e(r2_p)
local llike=e(ll)
test phintBW=phintWB
local p=r(p)
eststo m1: margins, dydx(phintBW phintWB highprob) post
estadd scalar p = `p'
estadd scalar r2p = `r2p'
estadd scalar llike = `llike'

eststo: logit ip_ phintBW phintWB highprob if blackhint==0, vce(cluster subject_id)
test phintBW=phintWB
local p=r(p)
local r2p=e(r2_p)
local llike=e(ll)

eststo m2: margins, dydx(phintBW phintWB highprob) post
estadd scalar p = `p'
estadd scalar r2p = `r2p'
estadd scalar llike = `llike'

eststo: logit ip_ phintBW phintWB highprob if blackhint==1, vce(cluster subject_id)
test phintBW=phintWB
local p=r(p)
local r2p=e(r2_p)
local llike=e(ll)
eststo m3: margins, dydx(phintBW phintWB highprob) post
estadd scalar p = `p'
estadd scalar r2p = `r2p'
estadd scalar llike = `llike'

eststo: logit ip_ phintBW phintWB highprob high_phintBW high_phintWB fef i.subject_id, vce(cluster subject_id)
test phintBW=phintWB
local p=r(p)
local r2p=e(r2_p)
local llike=e(ll)
eststo m4: margins, dydx(phintBW phintWB highprob high_phintBW high_phintWB fef) post
estadd scalar p = `p'
estadd scalar r2p = `r2p'
estadd scalar llike = `llike'

eststo: logit ip_ phintBW phintWB highprob high_phintBW high_phintWB fef i.subject_id if blackhint==0, vce(cluster subject_id)
test phintBW=phintWB
local p=r(p)
local r2p=e(r2_p)
local llike=e(ll)
eststo m5: margins, dydx(phintBW phintWB highprob high_phintBW high_phintWB fef) post
estadd scalar p = `p'
estadd scalar r2p = `r2p'
estadd scalar llike = `llike'

eststo: logit ip_ phintBW phintWB highprob high_phintBW high_phintWB fef i.subject_id if blackhint==1, vce(cluster subject_id)
test phintBW=phintWB
local p=r(p)
local r2p=e(r2_p)
local llike=e(ll)
eststo m6: margins, dydx(phintBW phintWB highprob high_phintBW high_phintWB fef) post
estadd scalar p = `p'
estadd scalar r2p = `r2p'
estadd scalar llike = `llike'

esttab m1 m2 m3 m4 m5 m6 using "./Tables/table_ip0.tex", b(%9.3g) t(%9.1f) ar2(%9.2f) stats(p N r2p llike, labels("P(FP rate $\neq$ FN rate)" "N" "Pseudo R-squared" "Log-likelihood")) indicate(Subject FE = fef) label addnotes("Errors are clustered by subject, average marginal treatment effects") mtitles("All" "S=White" "S=Black" "All" "S=White" "W=Black") title(Informed protection response: logistic regression) star("*" 0.10 "**" 0.05 "***" 0.01) nobaselevels compress nogaps replace



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
esttab m1 m2 m3 m4 using "./Tables/table_ip5.tex", b(%9.3g) t(%9.1f) ar2(%9.2f) label addnotes("Reporting average marginal effects, subject FE, errors are clustered by subject." "With flexible controls of posterior probability" ///
 ) mtitles("" "" "" "") title(Informed Protection Response: logit with flexible control for posteriors) star("*" 0.10 "**" 0.05 "***" 0.01) nobaselevels compress nogaps replace

 
*Saving the file to do heterogeneity analysis with respect to IP strategies:
save "./Temp/prep_heter.dta", replace


use "./Temp/prep_heter.dta", replace
*Test average response vs posterior probability by 




**Testing heterogeneity in informed protection response by signal color (black/white) and prior:
set more off
eststo clear
eststo: probit ip_ phintBW phintWB if blackhint==0&plevel==100, vce(cluster subject_id)
eststo: probit ip_ phintBW phintWB if blackhint==0&plevel==200, vce(cluster subject_id)
eststo: probit ip_ phintBW phintWB if blackhint==0&plevel==300, vce(cluster subject_id)
eststo: probit ip_ phintBW phintWB if blackhint==0&plevel==500, vce(cluster subject_id)

eststo: probit ip_ phintBW phintWB  if blackhint==1&plevel==100, vce(cluster subject_id)
eststo: probit ip_ phintBW phintWB  if blackhint==1&plevel==200, vce(cluster subject_id)
eststo: probit ip_ phintBW phintWB  if blackhint==1&plevel==300, vce(cluster subject_id)
eststo: probit ip_ phintBW phintWB  if blackhint==1&plevel==500, vce(cluster subject_id)

esttab using "./Tables/table_ip_het.tex", b(%9.3g) t(%9.1f) ar2(%9.2f) label addnotes("First four for white signal, the rest - black" ///
  "Errors are clustered by subject") mtitles("0.1" "0.2" "0.3" "0.5" "0.1" "0.2" "0.3" "0.5") title(Informed protection by prior) star("*" 0.10 "**" 0.05 "***" 0.01) nobaselevels compress nogaps replace

 
*Is it just an error in beliefs?
*Doing the IP regression and controlling for beliefs:
gen bes=be_^2
mkspline bespline1 0.2 bespline2 0.4 bespline3 0.6 bespline4 0.8 bespline5 = be_

  
  
**Robustness: flexible control both for beliefs and posteriors:
eststo clear
eststo: logit ip_ post_prob0* bespline* phintBW phintWB, vce(cluster subject_id)
eststo m1: margins, dydx(phintBW phintWB) post
eststo: logit ip_ post_prob0* bespline* phintBW phintWB i.subject_id, vce(cluster subject_id)
eststo m2: margins, dydx(phintBW phintWB) post
eststo: logit ip_ post_prob0* bespline* highprob phintBW phintWB high_phintBW high_phintWB i.subject_id, vce(cluster subject_id)
eststo m3: margins, dydx(highprob phintBW phintWB high_phintBW high_phintWB) post
eststo: logit ip_ post_prob0* bespline* blackhint phintBW phintWB black_phintBW black_phintWB i.subject_id, vce(cluster subject_id)
eststo m4: margins, dydx(blackhint phintBW phintWB black_phintBW black_phintWB) post
eststo: logit ip_ post_prob0* bespline* phintBW phintWB if blackhint==0, vce(cluster subject_id)
eststo m5: margins, dydx(phintBW phintWB) post
eststo: logit ip_ post_prob0* bespline* phintBW phintWB if blackhint==1, vce(cluster subject_id)
eststo m6: margins, dydx(phintBW phintWB) post
esttab m1 m2 m3 m4 m5 m6 using "./Tables/table_ip8_be.tex", b(%9.3g) t(%9.1f) ar2(%9.2f) label addnotes("With flexible controls of posterior probability and beliefs" ///
  "Errors are clustered by subject, average marginal treatment effects") mtitles("" "FE" "" "" "S=White" "S=Black") title(Informed Protection Response: flexible control for posteriors and beliefs) star("*" 0.10 "**" 0.05 "***" 0.01) nobaselevels compress nogaps replace

  
  
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

esttab m1 m2 m3 m4 using "./Tables/table_ip_flex.tex", b(%9.3g) t(%9.1f) label addnotes("With flexible controls of posterior probability and beliefs" ///
  "Subject FE, errors are clustered by subject, average marginal treatment effects") stats(N r2p llike, labels("N" "Pseudo R-squared" "Log-likelihood")) mtitles("Posterior only" "Posterior only" "Both" "Both") title(Informed Protection Response: flexible control for posteriors and beliefs) star("*" 0.10 "**" 0.05 "***" 0.01) nobaselevels compress nogaps replace


  
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

*Final robustness: semiparametric control for posteriors
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

merge m:1 subject_id using "./Temp/ipclasses.dta" //risk aversion, blind prot choices and demographic vars
drop _merge

gen class_alt=2-class
tab class class_alt
label var class_alt "IP strategy class"
label define clsnames 0 "Bayesians" 1 "Simpletons"
gen false_prob=phintWB+phintBW
label var false_prob "Prop. of lying gremlins"
label values class_alt clsnames



eststo clear
eststo: reg bel_err phintBW phintWB i.subject_id, vce(cluster subject_id)
eststo: reg bel_err phintBW phintWB i.subject_id if blackhint==0, vce(cluster subject_id)
eststo: reg bel_err phintBW phintWB i.subject_id if blackhint==1, vce(cluster subject_id)
esttab using "./Tables/table_be_err.tex", b(%9.3g) se(%9.1f) ar2(%9.2f) label addnotes(Dep. variable: reported belief - posterior probability) title(Belief Elicitation: When Mistakes Happen) mtitles("All" "S=White" "S=Black") star("*" 0.10 "**" 0.05 "***" 0.01) indicate(Subject FE = *.subject_id) nobaselevels compress nogaps replace


*Testing FP/FN confusion:
eststo clear
eststo: reg be_ phintBW phintWB i.plevel i.subject_id, vce(cluster subject_id)
eststo: reg be_ phintBW phintWB i.plevel if blackhint==0, vce(cluster subject_id)
eststo: reg be_ phintBW phintWB i.plevel if blackhint==1, vce(cluster subject_id)
esttab using "./Tables/table_be_err1.tex", b(%9.3g) se(%9.1f) ar2(%9.2f) label addnotes(Dep. variable: reported belief - posterior probability) title(Belief Elicitation: When Mistakes Happen) mtitles("All" "S=White" "S=Black") star("*" 0.10 "**" 0.05 "***" 0.01) indicate(Subject FE = *.subject_id) nobaselevels compress nogaps replace







eststo clear
eststo: reg bel_err i.class_alt##c.phintWB i.class_alt##c.phintBW i.subject_id, vce(cluster subject_id)
eststo: reg bel_err i.class_alt##c.phintWB i.class_alt##c.phintBW i.subject_id if blackhint==0, vce(cluster subject_id)
eststo: reg bel_err i.class_alt##c.phintWB i.class_alt##c.phintBW i.subject_id if blackhint==1, vce(cluster subject_id)
esttab using "./Tables/table_be_errx.tex", b(%9.3g) se(%9.1f) ar2(%9.2f) label addnotes(Dep. variable: reported belief - posterior probability) title(Belief Elicitation: When Mistakes Happen) mtitles("All" "S=White" "S=Black") star("*" 0.10 "**" 0.05 "***" 0.01) indicate(Subject FE = *.subject_id) nobaselevels compress nogaps replace



*Testing the accuracy of reported beliefs***
eststo clear
eststo: reg bel_err phintWB phintBW, vce(cluster subject_id)
eststo: reg bel_err i.plevel phintWB phintBW, vce(cluster subject_id)
eststo: reg bel_err i.goodquiz##c.phintWB i.goodquiz##c.phintBW, vce(cluster subject_id)
eststo: reg bel_err i.goodquiz##i.plevel i.goodquiz##c.phintWB i.goodquiz##c.phintBW, vce(cluster subject_id)
eststo: reg bel_err i.stat_educ##c.phintWB i.stat_educ##c.phintBW, vce(cluster subject_id)
eststo: reg bel_err i.stat_educ##i.plevel i.stat_educ##c.phintWB i.stat_educ##c.phintBW, vce(cluster subject_id)
esttab using "./Tables/table_be2.tex", b(%9.3g) se(%9.1f) ar2(%9.2f) label title(Belief Elicitation: Discrepancy) mtitles("" "" "" "" "" "") star("*" 0.10 "**" 0.05 "***" 0.01) indicate(Prior prob dummies = *.plevel) nobaselevels compress nogaps replace


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
esttab using "./Tables/table_be3.tex", b(%9.3g) t(%9.1f) ar2(%9.2f) label drop(_cons) mtitles("OLS" "FE" "OLS" "FE" "OLS" "FE") star("*" 0.10 "**" 0.05 "***" 0.01) note("Decomposition works only for imperfect signals") nobaselevels compress nogaps replace

gen psignalblack=p*phintBB+(1-p)*phintBW
gen psignalwhite=1-psignalblack
gen pifsignblack=p*phintBB/psignalblack

local c=0.43
local d=0.25
gen pseudo_psign=phintBB^`c'*p^`d'+phintBW^`c'*(1-p)^`d'


**********************************************
****-- Main Regressions: WTP FOR INFORMATION --****
**********************************************
set more off
clear all
use "./Output/main_waves.dta", replace
xtset subject_id round

merge m:1 subject_id round using "./Temp/bel_accuracy.dta" //average belief accuracy by subject
drop _merge


merge m:1 participant_id p using "./Temp/bp_val.dta" //risk aversion, blind prot choices and demographic vars
drop if _merge==2
drop _merge

merge m:1 subject_id using "./Temp/ipclasses.dta" //risk aversion, blind prot choices and demographic vars
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

*drop p_wei false_pos_wei false_neg_wei cost_bp_wei



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
label value accur_bel2 accur_bel_l


*Histograms:
hist ip_val_diff, title("Distribution of expected costs discrepancies") xtitle("Discrepancy") fraction note("Difference between optimal and actual expected costs in the IP treatment (by round)") color(navy)
graph export "./Graphs/hist_costs_discr.png", width(1200) height(800) replace



hist wtp, title("Distribution of WTP") xtitle("USD") fraction color(navy)
graph export "./Graphs/hist_WTP.png", width(1200) height(800) replace


hist value, title("Distribution of signal's value for a risk-neutral subject") xtitle("USD") fraction note("Value=expected change in costs from BP to IP") color(navy)
graph export "./Graphs/hist_value.png", width(1200) height(800) replace


hist wtp_diff, title("Distribution of WTP discrepancies (WTP - Value)") xtitle("USD") fraction note("Difference between stated wtp and theoretical value for a risk-neutral subject (each choice=obs)") color(navy)
graph export "./Graphs/hist_WTP_discr1.png", width(1200) height(800) replace



*Plotting discrepancies by signal type:
hist wtp_diff if fp_env==0&fn_env==0, title("Distribution of WTP discrepancies (WTP - Value)") xtitle("USD") fraction note("Honest signals") color(navy)
graph export "./Graphs/hist_WTP_discr1h.png", width(1200) height(800) replace


hist wtp_diff if fp_env==1&fn_env==0, title("Distribution of WTP discrepancies (WTP - Value)") xtitle("USD") fraction note("False-positive only signals") color(navy)
graph export "./Graphs/hist_WTP_discr1fp.png", width(1200) height(800) replace

hist wtp_diff if fp_env==0&fn_env==1, title("Distribution of WTP discrepancies (WTP - Value)") xtitle("USD") fraction note("False-negative only signals") color(navy)
graph export "./Graphs/hist_WTP_discr1fn.png", width(1200) height(800) replace

hist wtp_diff if fp_env==1&fn_env==1, title("Distribution of WTP discrepancies (WTP - Value)") xtitle("USD") fraction note("Positive false-positive and false-negative rates") color(navy)
graph export "./Graphs/hist_WTP_discr1fpfn.png", width(1200) height(800) replace



sort subject_id
gen wtp_diff_abs=abs(wtp_diff)
by subject_id: egen totwtp_diff=sum(wtp_diff_abs)
replace totwtp_diff=(1/6)*totwtp_diff

hist totwtp_diff, title("Distribution of WTP discrepancies (WTP - Value)") xtitle("USD") fraction note("Average absolute deviation between stated wtp and theoretical value for a risk-neutral subject, by subject") color(navy)
graph export "./Graphs/hist_WTP_discr2.png", width(1200) height(800) replace


gen highprob=p>0.201
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
label define risk_prefl 0 "Risk-neutral" 1 "Risk-loving" 2 "Risk-averse" 3 "No risk av. measure"
label value risk_pref risk_prefl

gen accur_bel1=1-accur_bel
label var accur_bel1 "Beliefs accuracy"
label def accur_bel1l 0 "Accur. beliefs" 1 "Inaccurate beliefs"
label value accur_bel1 accur_bel1l

replace accur_bel2=1-round(accur_bel2)
label value accur_bel2 accur_bel1l

*Self-reported strategies
*gen strategy_short=strategy_ip
*replace strategy_short=3 if strategy_short>3
label var strategy_short "Strategy"
label define strategy_short_l 1 "Rational" 2 "Seek honest" 3 "Other"
label values strategy_short strategy_short_l


label var be_change "Belief change"
label var confid "Certainty"



**Baseline WTP regs: tobit
eststo clear
eststo: tobit wtp false_neg false_pos, ll(0) ul(5)
eststo: tobit wtp false_neg false_pos cost_bp, ll(0) ul(5)
eststo: tobit wtp i.risk_pref##c.false_pos i.risk_pref##c.false_neg, ll(0) ul(5)
eststo: tobit wtp i.accur_bel2##c.false_pos i.accur_bel2##c.false_neg, ll(0) ul(5)
eststo: tobit wtp i.plevel##c.false_pos i.plevel##c.false_neg, ll(0) ul(5)
esttab using "./Tables/table_wtp_01tob.tex", b(%9.3g) se(%9.1f) ar2(%9.2f) label title(WTP for Information (tobit)) mtitles("" "" "" "" "") star("*" 0.10 "**" 0.05 "***" 0.01) nobaselevels compress nogaps replace


**Baseline WTP regs: tobit
gen class_alt=2-class

tab class class_alt
label var class_alt "IP strategy class"
label define clsnames 0 "Bayesians" 1 "Simpletons"
label values class_alt clsnames


eststo clear
eststo: tobit wtp false_neg false_pos, ll(0) ul(5)
eststo: tobit wtp i.class_alt##c.false_neg i.class_alt##c.false_pos, ll(0) ul(5)
eststo: tobit wtp false_neg false_pos, ll(0) ul(5)
eststo: tobit wtp i.class_alt##c.false_neg i.class_alt##c.false_pos, ll(0) ul(5)
esttab using "./Tables/table_wtp_classes.tex", b(%9.3g) se(%9.1f) ar2(%9.2f) label title(WTP for Information: heterogeneity by IP class) mtitles("p$<$0.3" "p$<$0.3" "All" "All") star("*" 0.10 "**" 0.05 "***" 0.01) nobaselevels compress nogaps replace


gen resp_freqBW=ip_b*freqBW
gen resp_freqBB=ip_b*freqBB
gen resp_freqWB=ip_w*freqWB

gen resp_freqBW1=(1-ip_b)*freqBW
gen resp_freqBB1=(1-ip_b)*freqBB
gen resp_freqWB1=(1-ip_w)*freqWB

reg wtp freqBW freqBB freqWB

reg wtp resp_freqBW resp_freqBB resp_freqWB

reg wtp resp_freqBW resp_freqBB resp_freqWB resp_freqBW1 resp_freqBB1 resp_freqWB1



**Baseline WTP regs 2: tobit
eststo clear
eststo: tobit wtp false_neg false_pos, ll(0) ul(5)
eststo: tobit wtp false_neg false_pos if plevel==100, ll(0) ul(5)
eststo: tobit wtp false_neg false_pos if plevel==200, ll(0) ul(5)
eststo: tobit wtp false_neg false_pos cost_bp, ll(0) ul(5)
eststo: tobit wtp false_neg false_pos cost_bp be_change, ll(0) ul(5)
eststo: tobit wtp false_neg false_pos cost_bp confid, ll(0) ul(5)
esttab using "./Tables/table_wtp_02tob.tex", b(%9.3g) se(%9.1f) ar2(%9.2f) label title(WTP for Information (tobit)) mtitles("All" "p=0.1" "p=0.2" "All" "All" "All") star("*" 0.10 "**" 0.05 "***" 0.01) nobaselevels compress nogaps replace


**Strategy effects
eststo clear
eststo: reg wtp_diff i.strategy_short false_neg false_pos
eststo: reg wtp_diff i.strategy_short##c.false_neg i.strategy_short##c.false_pos
eststo: reg wtp false_neg false_pos if plevel==100
eststo: reg wtp_diff i.strategy_short##c.false_neg i.strategy_short##c.false_pos if plevel==100
eststo: reg wtp false_neg false_pos if plevel==200
eststo: reg wtp_diff i.strategy_short##c.false_neg i.strategy_short##c.false_pos if plevel==200
esttab using "./Tables/table_wtp_strat.tex", b(%9.3g) se(%9.1f) ar2(%9.2f) label title(WTP minus Value of Information, connection to self-reported protection strategy) mtitles("All" "p=0.1" "p=0.2" "All" "All" "All") star("*" 0.10 "**" 0.05 "***" 0.01) nobaselevels compress nogaps replace



*Analyzing the effects of risk aversion on sensitivity of wtp to priors
eststo clear
eststo: reg wtp_diff i.highprob##c.false_neg i.highprob##c.false_pos, vce(cluster subject_id)
eststo: reg wtp_diff i.highprob##c.false_neg i.risk_pref#i.highprob#c.false_neg i.highprob##c.false_pos i.risk_pref#i.highprob#c.false_pos, vce(cluster subject_id)
eststo: reg wtp_diff i.risk_pref##i.highprob##c.false_neg i.risk_pref##i.highprob##c.false_pos, vce(cluster subject_id)
*eststo: reg wtp_diff i.highprob##c.false_neg c.totprot#i.highprob#c.false_neg i.highprob##c.false_pos c.totprot#i.highprob#c.false_pos, vce(cluster subject_id)

eststo: reghdfe wtp_diff i.highprob##c.false_neg i.risk_pref#i.highprob#c.false_neg i.highprob##c.false_pos i.risk_pref#i.highprob#c.false_pos, abs(subject_id) vce(cluster subject_id)
eststo: reghdfe wtp_diff i.risk_pref##i.highprob##c.false_neg i.risk_pref##i.highprob##c.false_pos, abs(subject_id) vce(cluster subject_id)
*eststo: reghdfe wtp_diff i.highprob##c.false_neg c.totprot#i.highprob#c.false_neg i.highprob##c.false_pos c.totprot#i.highprob#c.false_pos, abs(subject_id) vce(cluster subject_id)

esttab using "./Tables/wtp_het_risk.tex", b(%9.3g) se(%9.1f) ar2(%9.2f) label drop(_cons *.risk_pref#*.highprob *.risk_pref#c.false_pos *.risk_pref#c.false_neg) indicate("Full risk pref interactions=*.risk_pref") title(WTP minus Value of Information, risk aversion and sensitivity to FP and FN costs) mtitles("" "" "" "FE" "FE") star("*" 0.10 "**" 0.05 "***" 0.01) nobaselevels compress nogaps replace
esttab using "./Tables/wtp_het_risk_pres.tex", b(%9.3g) ar2(%9.2f) not label drop(_cons *.risk_pref#*.highprob *.risk_pref#*.highprob  *.risk_pref#c.false_pos *.risk_pref#c.false_neg) indicate("Full risk pref interactions=*.risk_pref") mtitles("" "" "" "FE" "FE") star("*" 0.10 "**" 0.05 "***" 0.01) nobaselevels compress nogaps replace



sort subject_id plevel
reghdfe wtp phintWB phintBW if plevel==100, abs(subject_id) vce(cluster subject_id)


*Analyzing changes in WTP across prior by subject to find if heterogeneous risk preferences can explain the pattern
save  "./Temp/wtp_discrepancy0.dta", replace

use "./Temp/wtp_discrepancy0.dta", replace
tab plevel phintBW
sort subject_id plevel round
keep subject_id plevel round wtp phintWB phintBW plevel value
reshape wide wtp phintWB phintBW plevel value, i(subject_id) j(round)

gen delta_b11=wtp1-wtp2
gen delta_b12=wtp4-wtp5

gen delta_b21=wtp1-wtp3
gen delta_b22=wtp4-wtp6

gen delta_v11=value1-value2
gen delta_v12=value4-value5

gen delta_v21=value1-value3
gen delta_v22=value4-value6

sum delta_b11 delta_b12
sum delta_b21 delta_b22

ttest delta_b11==delta_b12
ttest delta_b21==delta_b22

gen fp_env2=phintBW2>0
gen fn_env2=phintWB2>0

gen fp_env3=phintBW3>0
gen fn_env3=phintWB3>0
tab fp_env2 fn_env2

ttest delta_b11==delta_b12 if fn_env2==1
ttest delta_v11==delta_v12 if fn_env2==1

ttest delta_b11==delta_b12 if fp_env2==1
ttest delta_v11==delta_v12 if fp_env2==1

ttest delta_b21==delta_b22 if fn_env3==1
ttest delta_v21==delta_v22 if fn_env3==1

ttest delta_b11==delta_b12 if fp_env2==1
ttest delta_v11==delta_v12 if fp_env2==1


stop

gen delta_b11x=(wtp[_n-1]-wtp)-1000*(round~=2)
gen delta_b12x=(wtp[_n-1]-wtp)-1000*(round~=5)
gen delta_b21x=(wtp[_n-2]-wtp)-1000*(round~=3)
gen delta_b22x=(wtp[_n-2]-wtp)-1000*(round~=6)


by subject_id: egen delta_b11=max(delta_b11x)
by subject_id: egen delta_b21=max(delta_b21x)
by subject_id: egen delta_b12=max(delta_b12x)
by subject_id: egen delta_b22=max(delta_b22x)

gen delta_b11x=(wtp[_n-1]-wtp)-1000*(round~=2)
gen delta_b12x=(wtp[_n-1]-wtp)-1000*(round~=5)
gen delta_b21x=(wtp[_n-2]-wtp)-1000*(round~=3)
gen delta_b22x=(wtp[_n-2]-wtp)-1000*(round~=6)


by subject_id: egen delta_b11=max(delta_b11x)
by subject_id: egen delta_b21=max(delta_b21x)
by subject_id: egen delta_b12=max(delta_b12x)
by subject_id: egen delta_b22=max(delta_b22x)




sum delta_b1 delta_b2 delta_v1 delta_v2
order subject_id round plevel phintWB phintBW wtp delta_b1 delta_b2
drop delta_b1x delta_b2x delta_v1x delta_v2x
bro

use "./Temp/wtp_discrepancy0.dta", replace

gen fp_cost=phintBW*protectioncost
gen fn_cost=phintWB*loss


gen phintWBcat=round(100*phintWB)
gen phintBWcat=round(100*phintBW)

gen false_cost=fp_cost+fn_cost

reg wtp i.plevel false_cost
reg wtp i.plevel fp_cost fn_cost
reg wtp i.plevel i.risk_pref##c.fp_cost i.risk_pref##c.fn_cost


gen delta_p=p-0.275
gen false_pos2=delta_p*fp_cost
gen false_neg2=delta_p*fn_cost

label var fp_cost "FP rate"
label var fn_cost "FN rate"
label var false_pos2 "Prior_change$\times\$FP rate"
label var false_neg2 "Prior_change$\times\$FN rate"


*Testing if subjects' sensitivity to FP/FN costs changes with priors:
eststo clear
*eststo: reg wtp i.plevel fp_cost fn_cost, vce(cluster subject_id)
*estadd scalar p = 1
*estadd scalar adj_r= e(r2_a)
eststo: reg wtp i.plevel fp_cost fn_cost false_pos2 false_neg2, vce(cluster subject_id)
estadd scalar adj_r= e(r2_a)
test false_pos2 false_neg2
local p=r(p)
estadd scalar p = `p'
eststo: reg wtp i.plevel fp_cost fn_cost i.plevel#c.fp_cost i.plevel#c.fn_cost, vce(cluster subject_id)
estadd scalar adj_r= e(r2_a)
test 200.plevel#c.fp_cost 200.plevel#c.fn_cost 300.plevel#c.fp_cost 300.plevel#c.fn_cost 500.plevel#c.fp_cost 500.plevel#c.fn_cost
local p=r(p)
estadd scalar p = `p'
*eststo: reghdfe wtp i.plevel fp_cost fn_cost, absorb(subject_id) vce(cluster subject_id)
*estadd scalar adj_r= e(r2_a)
*estadd scalar p = 1
eststo: reghdfe wtp i.plevel fp_cost fn_cost false_pos2 false_neg2, absorb(subject_id) vce(cluster subject_id)
estadd scalar adj_r= e(r2_a)
test false_pos2 false_neg2
local p=r(p)
estadd scalar p = `p'
eststo: reghdfe wtp i.plevel fp_cost fn_cost i.plevel#c.fp_cost i.plevel#c.fn_cost, absorb(subject_id) vce(cluster subject_id)
estadd scalar adj_r= e(r2_a)
test 200.plevel#c.fp_cost 200.plevel#c.fn_cost 300.plevel#c.fp_cost 300.plevel#c.fn_cost 500.plevel#c.fp_cost 500.plevel#c.fn_cost
local p=r(p)
estadd scalar p = `p'
esttab using "./Tables/wtp_het_test.tex", b(%9.3g) se(%9.1f) ar2(%9.2f) label drop(_cons fp_cost fn_cost *.plevel) stats(p adj_r N, labels("P(beta_X=0)" "Adjusted \(R^{2}\)" "Observations")) title(WTP for Information (OLS): analyze FP/FN sensitivity by prior) mtitles("" "" "" "FE" "FE" "FE") star("*" 0.10 "**" 0.05 "***" 0.01) nobaselevels compress nogaps replace
esttab using "./Tables/wtp_het_test_pres.tex", b(%9.3g) se(%9.1f) ar2(%9.2f) label drop(_cons fp_cost fn_cost *.plevel) stats(p adj_r N, labels("P(beta_X=0)" "Adjusted \(R^{2}\)" "Observations")) addnotes("Controlling for FP/FN rates, prior, constant omitted") mtitles("" "" "" "FE" "FE" "FE") star("*" 0.10 "**" 0.05 "***" 0.01) nobaselevels compress nogaps replace


*Testing if there is distinguishing btw FP and FN rates:
eststo clear
eststo: reg wtp i.plevel phintWB phintBW, vce(cluster subject_id)
test phintWB=phintBW
estadd scalar adj_r= e(r2_a)
estadd scalar p = r(p)
eststo: reg wtp i.plevel i.plevel#c.phintWB i.plevel#c.phintBW, vce(cluster subject_id)
test  (100.plevel#c.phintWB=100.plevel#c.phintBW) (200.plevel#c.phintWB=200.plevel#c.phintBW)  (300.plevel#c.phintWB=300.plevel#c.phintBW)   (500.plevel#c.phintWB=500.plevel#c.phintBW)
estadd scalar adj_r= e(r2_a)
estadd scalar p = r(p)
eststo: reghdfe wtp i.plevel phintWB phintBW,  absorb(subject_id) vce(cluster subject_id)
test phintWB=phintBW
estadd scalar adj_r= e(r2_a)
estadd scalar p = r(p)
eststo: reghdfe wtp i.plevel i.plevel#c.phintWB i.plevel#c.phintBW,  absorb(subject_id) vce(cluster subject_id)
test  (100.plevel#c.phintWB=100.plevel#c.phintBW) (200.plevel#c.phintWB=200.plevel#c.phintBW)  (300.plevel#c.phintWB=300.plevel#c.phintBW)   (500.plevel#c.phintWB=500.plevel#c.phintBW)
estadd scalar adj_r= e(r2_a)
estadd scalar p = r(p)
esttab using "./Tables/wtp_fpvsfn_test.tex", b(%9.3g) se(%9.1f) ar2(%9.2f) label drop(_cons *.plevel) stats(p adj_r N, labels("P(FP rate=FN rate)" "Adjusted \(R^{2}\)" "Observations")) title(WTP for Information (OLS): comparing sensitivity to FP/FN rates) mtitles("" "" "FE" "FE") star("*" 0.10 "**" 0.05 "***" 0.01) nobaselevels compress nogaps replace
esttab using "./Tables/wtp_fpvsfn_test_pres.tex", b(%9.3g) se(%9.1f) ar2(%9.2f) label drop(_cons *.plevel) stats(p adj_r N, labels("P(FP rate=FN rate)" "Adjusted \(R^{2}\)" "Observations")) mtitles("" "" "FE" "FE") star("*" 0.10 "**" 0.05 "***" 0.01) nobaselevels compress nogaps replace





gen false_prob=phintWB+phintBW


gen false_tot=false_pos+false_neg


**Baseline WTP difference: risk-aversion and belief accuracy:
eststo clear
eststo: reg wtp_diff false_pos false_neg, vce(cluster subject_id)
eststo: reghdfe wtp_diff false_pos false_neg, abs(subject_id) vce(cluster subject_id)
eststo: reghdfe wtp_diff i.risk_pref##c.false_pos i.risk_pref##c.false_neg, abs(subject_id) vce(cluster subject_id)
eststo: reghdfe wtp_diff i.accur_bel2##c.false_pos i.accur_bel2##c.false_neg, abs(subject_id) vce(cluster subject_id)
eststo: reghdfe wtp_diff i.plevel##c.false_pos i.plevel##c.false_neg, abs(subject_id) vce(cluster subject_id)
*eststo: reghdfe wtp_diff i.accur_bel2#i.risk_pref##c.false_pos i.accur_bel2#i.risk_pref##c.false_neg if plevel<300, abs(subject_id) vce(cluster subject_id)
esttab using "./Tables/table_wtpdiff_01rfull_ag.tex", b(%9.3g) se(%9.1f) ar2(%9.2f) label title(WTP minus Value of Information (OLS)) mtitles("" "" "" "" "") star("*" 0.10 "**" 0.05 "***" 0.01) nobaselevels compress nogaps replace
esttab using "./Tables/table_wtpdiff_01rfull_ag_pres.tex", b(%9.3g) ar2(%9.2f) not indicate("Risk pref dummies = *.risk_pref" "Prior dummies = *.plevel" "Belief accuracy = *.accur_bel2" ) label mtitles("" "FE" "FE" "FE" "FE") star("*" 0.10 "**" 0.05 "***" 0.01) nobaselevels compress nogaps replace
esttab using "./Tables/table_wtpdiff_01r_ag.tex", b(%9.3g) se(%9.1f) ar2(%9.2f) label title(WTP minus Value of Information (OLS)) mtitles("" "" "" "" "") star("*" 0.10 "**" 0.05 "***" 0.01) noomit nobaselevels compress nogaps replace





*Demographic variables:
eststo clear
*eststo: reg wtp_diff false_pos false_neg, vce(cluster subject_id)
eststo: reg wtp_diff false_pos false_neg i. sex i.sex#c.false_pos i.sex#c.false_neg, vce(cluster subject_id)
eststo: reg wtp_diff i.sex##i.plevel i.sex##c.false_pos i.sex##c.false_neg, vce(cluster subject_id)
eststo: reg wtp_diff i.stat_educ##c.false_pos i.stat_educ##c.false_neg, vce(cluster subject_id)
eststo: reg wtp_diff i.stat_educ##i.plevel i.stat_educ##c.false_pos i.stat_educ##c.false_neg, vce(cluster subject_id)
eststo: reg wtp_diff i.old##c.false_pos i.old##c.false_neg, vce(cluster subject_id)
eststo: reg wtp_diff i.old##i.plevel i.old##c.false_pos i.old##c.false_neg, vce(cluster subject_id)
esttab using "./Tables/table_wtpdiff_02.tex", b(%9.3g) se(%9.1f) ar2(%9.2f) label indicate(Prior dummies=*.plevel) title(WTP minus Value of Information: demographic determinants) mtitles("" "" "" "" "" "" "" "" "") star("*" 0.10 "**" 0.05 "***" 0.01) nobaselevels compress nogaps replace


*Prior heterogeneity and instructions understanding:
eststo clear
eststo: reg wtp_diff i.plevel i.highprob##c.false_neg i.goodquiz#i.highprob#c.false_neg i.highprob##c.false_pos i.goodquiz#i.highprob#c.false_pos, vce(cluster subject_id)
eststo: reg wtp_diff i.plevel i.goodquiz##i.highprob##c.false_neg i.goodquiz##i.highprob##c.false_pos, vce(cluster subject_id)
*eststo: reg wtp_diff i.highprob##c.false_neg c.totprot#i.highprob#c.false_neg i.highprob##c.false_pos c.totprot#i.highprob#c.false_pos, vce(cluster subject_id)

eststo: reghdfe wtp_diff i.plevel i.highprob##c.false_neg i.goodquiz#i.highprob#c.false_neg i.highprob##c.false_pos i.goodquiz#i.highprob#c.false_pos, abs(subject_id) vce(cluster subject_id)
eststo: reghdfe wtp_diff i.plevel i.goodquiz##i.highprob##c.false_neg i.goodquiz##i.highprob##c.false_pos, abs(subject_id) vce(cluster subject_id)

esttab using "./Tables/table_wtpdiff_instr.tex", b(%9.3g) se(%9.1f) ar2(%9.2f) label indicate(Prior dummies=*.plevel) title(WTP minus Value of Information: instructions understanding) mtitles("" "" "" "" "" "" "" "" "") star("*" 0.10 "**" 0.05 "***" 0.01) nobaselevels compress nogaps replace



eststo clear
eststo: reg wtp_diff i.highprob##c.false_neg i.highprob##c.false_pos if goodquiz==0, vce(cluster subject_id)
eststo: reg wtp_diff i.highprob##c.false_neg i.highprob##c.false_pos if goodquiz==1, vce(cluster subject_id)

eststo: reg wtp_diff i.plevel i.highprob##c.false_neg i.highprob##c.false_pos if goodquiz==0, vce(cluster subject_id)
eststo: reg wtp_diff i.plevel i.highprob##c.false_neg i.highprob##c.false_pos if goodquiz==1, vce(cluster subject_id)
*eststo: reg wtp_diff i.highprob##c.false_neg c.totprot#i.highprob#c.false_neg i.highprob##c.false_pos c.totprot#i.highprob#c.false_pos, vce(cluster subject_id)
esttab using "./Tables/table_wtpdiff_instr.tex", b(%9.3g) se(%9.1f) ar2(%9.2f) label indicate(Prior dummies=*.plevel) title(WTP minus Value of Information: instructions understanding) mtitles("Bad" "Good" "Bad" "Good") star("*" 0.10 "**" 0.05 "***" 0.01) nobaselevels compress nogaps replace






*By prior prob for FP and FN costs:
eststo clear
eststo: reghdfe wtp_diff false_pos false_neg if plevel==100, abs(subject_id) vce(cluster subject_id)
eststo: reghdfe wtp_diff false_pos false_neg if plevel==200, abs(subject_id) vce(cluster subject_id)
eststo: reghdfe wtp_diff false_pos false_neg if plevel==300, abs(subject_id) vce(cluster subject_id)
eststo: reghdfe wtp_diff false_pos false_neg if plevel==500, abs(subject_id) vce(cluster subject_id)
esttab using "./Tables/table_wtpdiff_03.tex", b(%9.3g) se(%9.1f) ar2(%9.2f) label title(WTP - Value of Information, by prior) addnotes("Subject fixed effects are included.") mtitles("0.1" "0.2" "0.3" "0.5") star("*" 0.10 "**" 0.05 "***" 0.01) nobaselevels compress nogaps replace


*By prior prob for FP and FN rates:
eststo clear
eststo: tobit wtp phintBW phintWB if plevel==100, ll(0) ul(5)
test phintBW=phintWB
local p=r(p)
estadd scalar p = `p'
estadd scalar adj_r= e(r2_a)
eststo: tobit wtp phintBW phintWB if plevel==200, ll(0) ul(5)
test phintBW=phintWB
local p=r(p)
estadd scalar p = `p'
estadd scalar adj_r= e(r2_a)
eststo: tobit wtp phintBW phintWB if plevel==300, ll(0) ul(5)
test phintBW=phintWB
local p=r(p)
estadd scalar p = `p'
estadd scalar adj_r= e(r2_a)
eststo: tobit wtp phintBW phintWB if plevel==500, ll(0) ul(5)
test phintBW=phintWB
local p=r(p)
estadd scalar p = `p'
estadd scalar adj_r= e(r2_a)
esttab using "./Tables/table_wtpdiff_04tob.tex", b(%9.3g) se(%9.1f) ar2(%9.2f) label stats(p adj_r N, labels("P(FP rate=FN rate)" "Adjusted \(R^{2}\)" "Observations")) title(WTP for Information, by prior (tobit)) mtitles("0.1" "0.2" "0.3" "0.5") star("*" 0.10 "**" 0.05 "***" 0.01) nobaselevels compress nogaps replace
esttab using "./Tables/table_wtpdiff_04tob_pres.tex", b(%9.3g) se(%9.1f) ar2(%9.2f) drop(_cons) stats(p adj_r N, labels("P(FP rate=FN rate)" "Adjusted \(R^{2}\)" "Observations")) eqlabels(none) label mtitles("0.1" "0.2" "0.3" "0.5") addnotes("Tobit regression of WTP, constant omitted") star("*" 0.10 "**" 0.05 "***" 0.01) nobaselevels compress nogaps replace




save  "./Temp/wtp_discrepancy0.dta", replace

use "./Temp/wtp_discrepancy0.dta", replace
drop highprob
gen highprob=p>0.201


**Calculate average WTP discrepancy by signal type
tempname p1
postfile `p1' false_pos false_neg wtp_diff ptest using "./Temp/wtp_by_environment.dta", replace


forvalues i=0/1{
  forvalues j=0/1{
	    ttest wtp_diff==0 if fp_env==`i'&fn_env==`j', level(95)
		local ptest=r(p)
		local wtp_diff=r(mu_1)
		post `p1' (`i') (`j') (`wtp_diff') (`ptest')
    }
  }

postclose `p1'


**Calculate average WTP discrepancy by signal type
tempname p2
postfile `p2' highprob false_pos false_neg wtp_diff ptest using "./Temp/wtp_by_environment_det.dta", replace

forvalues k=0/1{
  forvalues i=0/1{
    forvalues j=0/1{
	    ttest wtp_diff==0 if highprob==`k'&fp_env==`i'&fn_env==`j', level(95)
		local ptest=r(p)
		local wtp_diff=r(mu_1)
		post `p2' (`k') (`i') (`j') (`wtp_diff') (`ptest')
    }
  }
}
postclose `p2'

use "./Temp/wtp_by_environment.dta", replace

tostring false_pos false_neg, replace
foreach var of varlist false_pos false_neg{
  replace `var'="No" if `var'=="0"
  replace `var'="Yes" if `var'=="1"
}
bro
format wtp_diff ptest %9.3f
listtex using "./Tables/bigpicture_wtp.tex", type rstyle(tabular) head("\begin{table}[H]\centering \caption{Average WTP discrepancy (WTP-Value) by Signal Type} \begin{tabular}{cccc} \hline \hline" `"\textbf{False-positive}&\textbf{False-negative}&\textbf{Mean WTP discrepancy}& \textbf{P($=0$)}\\ \hline"') foot("\hline \end{tabular} \end{table}") replace
listtex using "./Tables/bigpicture_wtp_pres.tex", type rstyle(tabular) head("\begin{table}[H]\centering \begin{tabular}{cccc} \hline \hline" `"\textbf{False-positive}&\textbf{False-negative}&\textbf{Mean WTP discrepancy}& \textbf{P($=0$)}\\ \hline"') foot("\hline \end{tabular} \end{table}") replace




use "./Temp/wtp_discrepancy0.dta", replace
*How many switch choices after changing prior probabilities:
gen wtp_p_change=wtp~=wtp[_n-3]
replace wtp_p_change=0 if round<4
bys subject_id: egen nchanges=sum(wtp_p_change)

gen wtp_s_change=wtp~=wtp[_n-1]
replace wtp_s_change=0 if inlist(round,1,2,4)
bys subject_id: egen nchanges_s=sum(wtp_s_change)
tab nchanges_s

**Removing subjects which do not update their WTP after increasing priors:
eststo clear
eststo: reghdfe wtp_diff phintBW phintWB if nchanges>0, abs(subject_id) vce(cluster subject_id)
eststo: reghdfe wtp_diff phintBW phintWB if plevel==100&nchanges>0, abs(subject_id) vce(cluster subject_id)
eststo: reghdfe wtp_diff phintBW phintWB if plevel==200&nchanges>0, abs(subject_id) vce(cluster subject_id)
eststo: reghdfe wtp_diff phintBW phintWB if plevel==300&nchanges>0, abs(subject_id) vce(cluster subject_id)
eststo: reghdfe wtp_diff phintBW phintWB if plevel==500&nchanges>0, abs(subject_id) vce(cluster subject_id)
esttab using "./Tables/table_wtpdiff_06s.tex", b(%9.3g) se(%9.1f) ar2(%9.2f) label title(WTP - Value of Information, by prior) mtitles("All" "0.1" "0.2" "0.3" "0.5")  addnotes("Only subjects who change their decisions across priors") star("*" 0.10 "**" 0.05 "***" 0.01) nobaselevels compress nogaps replace



*By prior prob with high/low prob: testing for the presentation order effect
*No variation in the presentation order. Comparing 0.1 vs 0.2 and 0.3 vs 0.5
gen highsession=seq>3
gen firstorder=round<4

label var highsession "Higher probabilities session"
label define highsessionl 0 "Low probabilities" 1 "Starts with p=0.2"
label value highsession highsessionl

label var firstorder "First prior in the sequence"
label define firstorderl 0 "Second prior" 1 "First prior"
label value firstorder firstorderl

tab highsession firstorder


eststo clear
eststo: reg wtp_diff phintBW phintWB if p<0.3, vce(cluster subject_id)
eststo: reg wtp_diff phintBW phintWB if p>0.25, vce(cluster subject_id)
eststo: reg wtp_diff i.highsession##c.phintBW i.highsession##c.phintWB if p<0.3, vce(cluster subject_id)
*eststo: reg wtp_diff phintBW phintWB if highsession==0, vce(cluster subject_id)
*eststo: reg wtp_diff phintBW phintWB if highsession==1, vce(cluster subject_id)
eststo: reg wtp_diff i.highsession##c.phintBW i.highsession##c.phintWB, vce(cluster subject_id)
eststo: reg wtp_diff i.firstorder##c.phintBW i.firstorder##c.phintWB, vce(cluster subject_id)
eststo: reg wtp_diff i.firstorder##c.phintBW i.firstorder##c.phintWB i.highsession#c.phintBW i.highsession#c.phintWB, vce(cluster subject_id)
esttab using "./Tables/table_wtpdiff_order.tex", b(%9.3g) se(%9.1f) ar2(%9.2f) label title(WTP - Value of Information, by prior with order effects) mtitles( "p=0.1,0.2" "p=0.3,0.5" "p=0.1,0.2" "" "" "") star("*" 0.10 "**" 0.05 "***" 0.01) nobaselevels compress nogaps replace






*Testing heterogeneity
reg wtp_diff i.plevel##c.false_pos i.plevel##c.false_neg, vce(cluster subject_id)
contrast plevel plevel#c.false_pos plevel#c.false_neg, overall


*By prior prob 2:
eststo clear
eststo: reghdfe  wtp_diff phintBW phintWB, abs(subject_id) vce(cluster subject_id)
eststo: reghdfe  wtp_diff phintBW phintWB if plevel==100, abs(subject_id) vce(cluster subject_id)
eststo: reghdfe  wtp_diff phintBW phintWB if plevel==200, abs(subject_id) vce(cluster subject_id)
eststo: reghdfe  wtp_diff phintBW phintWB if plevel==300, abs(subject_id) vce(cluster subject_id)
eststo: reghdfe  wtp_diff phintBW phintWB if plevel==500, abs(subject_id) vce(cluster subject_id)
esttab using "./Tables/table_wtpdiff_06.tex", b(%9.3g) se(%9.1f) ar2(%9.2f) label title(WTP - Value of Information, by prior) mtitles("All" "0.1" "0.2" "0.3" "0.5") star("*" 0.10 "**" 0.05 "***" 0.01) nobaselevels compress nogaps replace



*By prior prob:
eststo clear
eststo: tobit wtp highprob phintBW i.highprob#c.phintBW phintWB i.highprob#c.phintWB, ll(0) ul(5)
eststo: tobit wtp highprob phintBW i.highprob#c.phintBW phintWB i.highprob#c.phintWB if goodquiz==1, ll(0) ul(5)
eststo: tobit wtp highprob phintBW i.highprob#c.phintBW phintWB i.highprob#c.phintWB if stat_educ==1, ll(0) ul(5)
*eststo: tobit value highprob phintBW i.highprob#c.phintBW phintWB i.highprob#c.phintWB, ll(0) ul(5)
esttab using "./Tables/table_wtp_val0.tex", b(%9.3g) se(%9.1f) ar2(%9.2f) label title() mtitles("All" "Good quiz only" "Stat. classes") star("*" 0.10 "**" 0.05 "***" 0.01) nobaselevels compress nogaps replace



*Create variables for the next table:
gen ps2=p^2
label var ps2 "p squared"
label var freqBW "FP total prob."
label var freqWB "FN total prob."
gen freqBWs=freqBW^2
gen freqWBs=freqWB^2
label var freqBWs "FP total prob. sq."
label var freqWBs "FN total prob. sq."

*Interaction stuff:
eststo clear
eststo: tobit wtp p freqBW freqWB phintBW phintWB, ll(0) ul(5)
eststo: tobit wtp p freqBW freqWB phintBW phintWB ps2 freqBWs freqWBs, ll(0) ul(5)
eststo: tobit wtp p freqBW freqWB i.stat_educ##c.phintBW i.stat_educ##c.phintWB, ll(0) ul(5)
eststo: tobit wtp p freqBW freqWB i.stat_educ##c.phintBW i.stat_educ##c.phintWB ps2 freqBWs freqWBs, ll(0) ul(5)
*eststo: tobit value p freqBW freqWB phintBW phintWB, ll(0) ul(5)
*eststo: tobit value p freqBW freqWB phintBW phintWB ps2 freqBWs freqWBs, ll(0) ul(5)
esttab using "./Tables/table_wtp_val1.tex", b(%9.3g) se(%9.1f) ar2(%9.2f) label drop(p freqBW freqWB ps2 freqWBs) indicate(With squares=freqBWs) title() note("Controlling for priors and total probabilities of false-posiive and false-negative outcomes. Standard errors in parentheses.") mtitles("" "" "" "") star("*" 0.10 "**" 0.05 "***" 0.01) nobaselevels compress nogaps replace


*Paying positive amounts for signals not affecting their protection choices:
gen pay_notuse=(wtp>0)&(ip_b==ip_w)
gen use_signal=(ip_b>ip_w)
tab pay_notuse //many choices are inconsistent

//most subjects make at least one inconsistent choice here:
sort subject_id
by subject_id: egen totpay_notuse=sum(pay_notuse)

save "./Temp/wtpdat2.dta", replace


log close


