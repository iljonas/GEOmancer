@echo off
title GEO Assembler Runner
set Rbase="C:\Program Files\R"
for /F "tokens=*" %%i in ('dir %Rbase% /A:D /O:D /T:C /B') do set Rversion=%%i
set Rpath="%Rbase:"=%\%Rversion%\bin\Rscript.exe"

set scriptPath="C:\Users\%USERNAME%\Documents\Capstone Files\GEO-Antimicrobial-Adjunct-Project\GEO_Assembler.R"
set /p series="Enter Series ID: "
set /p master="Enter master file name (do not include file extension): "
set /p exclusion="Enter column exclusions (no extra spaces): "

echo Analyzing Series and Platform Files...
%Rpath% %scriptPath% %series% %master% %exclusion%

pause