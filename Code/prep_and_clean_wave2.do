**DATA PREPARATION AND CLEANING***
set more off
clear all


*!!put your working folder below:
*!!put your working folder below:*
cd <project_folder>

set seed 135

*Prepare the treatment characteristics to merge:
import delimited using "./Input/exp_treatments_wave2.csv", varnames(1) clear
*drop v1
rename ïsnames snames
rename snames treatn

save "./Temp/exp_treatments_wave2.dta", replace

import delimited using "./Input/treatment_sequences_wave2.csv", varnames(1) clear //I prepare this file separately in Excel, describes treatment order for each group of subjects/session type
*rename ïseq seq
reshape long r, i(seq) j(round) //so it is sequence - round structure to make merging with the results file easier
rename r treatn
*stop
merge m:1 treatn using "./Temp/exp_treatments_wave2.dta" //merging characteristics of each treatment (prior prob, info structure = N of each gremlin types)
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
drop if missing(treatn)
drop if missing(seq)

save "./Temp/exp_treatments_wave2.dta", replace


*First importing the Qualtrics data
*The results data has variable names incompatible with Stata
*Hence we import the first row of the results file first as data. The first row contains original variable names, which we process to an acceptable format and store in a string 
import delimited using "./Input/Data_Wave3.csv", rowrange(1:1) varnames(nonames) clear
local myvarlist
foreach v of varlist _all{
 display `v'
 replace `v'=strtrim(`v')
 replace `v'=subinstr(`v', " (in seconds)", "", .)
 replace `v'=subinstr(`v', " ", "", .)
 replace `v'=substr(`v',3,.)+"_"+substr(`v',1,1) if regexm(`v',"^[0-9]")
 replace `v'=strlower(`v')
 replace `v'=subinstr(`v', ".", "", .)
 replace `v' = substr(`v', 1, 32)
 local luck=`v'[1]
 local myvarlist `myvarlist' `luck'
}
*Then we import the data from the results file (ignoring first 3 rows with var names and definition). Then we rename these variables using the previously stored processed names 
import delimited using "./Input/Data_Wave3.csv", rowrange(4:106) varnames(nonames) clear
*stop
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

replace participant_id=participant_id+"w2"

rename q195_1_text gpa
replace gpa=. if gpa==-99


rename q210 acad_status
tostring acad_status, force replace
replace acad_status="Undergraduate student" if acad_status=="1"
replace acad_status="Graduate student" if acad_status=="2"
replace acad_status="Other" if acad_status=="3"

tab acad_status

save "./Temp/secondwave_wide.dta", replace

tab seq

*generating wave numbers as we will need them to disambiguate subject_id
gen wave=2

append using "./Temp/mainwaves_wide.dta", generate(_from)
tab seq
replace wave=1 if missing(wave)
gen student=((educ==3)|(educ==4))
save "./Temp/allwaves_wide.dta", replace
duplicates list participant_id
tab student if wave==1
tab student


*



*Exporting IP explanations for the IP (wave 2 as wave 1 is alredy coded):
keep if wave==2
keep participant_id q105
order participant_id q105
rename q105 ip_explanation
export delimited using "./Output/ip_explanations_wave2.csv", replace
drop if ip_explanation=="-99"
replace ip_explanation = subinstr(ip_explanation, "&", "\&", .)
replace ip_explanation = subinstr(ip_explanation, "%", "\%", .)
replace ip_explanation = subinstr(ip_explanation, "$", "\$", .)
replace ip_explanation="\item "+ip_explanation
keep ip_explanation

*drop if ip_explanation==-99
listtex using "./Output/wave2_explanations.txt", replace


*then exporting to tex to add to the paper:


use "./Temp/allwaves_wide.dta", replace
