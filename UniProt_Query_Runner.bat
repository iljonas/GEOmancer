@echo off
title UniProt Query Runner
set pyPath=C:\Users\%USERNAME%\AppData\Local\Programs\Python\Python36-32\python.exe

set Rbase="C:\Program Files\R"
for /F "tokens=*" %%i in ('dir %Rbase% /A:D /O:D /T:C /B') do set Rversion=%%i
set Rpath="%Rbase:"=%\%Rversion%\bin\Rscript.exe"

set sourcePath="C:\Users\%USERNAME%\Documents\Capstone Files\GEO-Antimicrobial-Adjunct-Project
set queryPath=%sourcePath%\UniProt_Query.py"
set transformPath=%sourcePath%\UniProt_Transformer.R"

set /p input="Enter source file name (do not include file extension): "

echo Running Expression IDs through UniProt
%pyPath% %queryPath% %input%

echo.
echo Combining UniProt output with project master file
%RPath% %transformPath% %input%

pause