**PRODUCES ALL THE TABLES AND GRAPHS FOR THE PROJECT**
*requires moremata, estout, reghdfe libraries
cd C:\Tornado_warnings\Experiment\Alerts_Experiment
do ".\Code\prep_and_clean.do"
do ".\Code\bp_analysis.do"
do ".\Code\summary_ip_be.do"
do ".\Code\base_regs_ip_be.do"
do ".\Code\main_regs_wtp.do"
do ".\Code\graph_sensitivities.do"
