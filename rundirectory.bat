:: producing all the outputs in the project
:: first deleting all the temporary files:
del /q ".\Temp\*"
del /q ".\Graphs\*"
del /q ".\Tables\*"
del /q ".\Output\*"
:: FOR /D %%p IN (".\Temp\*.*") DO rmdir "%%p" /s /q

:: running R code:
R CMD BATCH Code/analysis_treatments4.R
R CMD BATCH Code/firstpilotrun.R

:: Only Stata results (using the output from R scripts above)
Stata-64 /b do Code/pilot_analysis.do "./Temp/pilot_analysis.log"