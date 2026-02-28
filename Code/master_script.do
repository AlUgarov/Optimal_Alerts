**PRODUCES ALL THE TABLES AND GRAPHS FOR THE PROJECT**
*requires moremata, estout, ranktest, egenmore, tabout, listtex, aaplot, reghdfe, ivreg2 libraries
*cd <your folder project folder>
do ".\Code\prep_and_clean_wave1.do"
do ".\Code\prep_and_clean_wave2.do"
do ".\Code\bp_analysis_all.do"
do ".\Code\summary_ip_be_all_weighted.do"
do ".\Code\base_regs_ip_be_all_treatcl.do"
do ".\Code\main_regs_wtp_all_treatcl.do"
do ".\Code\graph_sensitivities_all.do"
do ".\Code\demogr_tables.do"
