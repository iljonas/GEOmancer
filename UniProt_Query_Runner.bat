@echo off
title UniProt Query Runner
set pyPath=C:\Users\%USERNAME%\AppData\Local\Programs\Python\Python36-32\python.exe

set scriptPath="C:\Users\%USERNAME%\Documents\Capstone Files\GEO-Antimicrobial-Adjunct-Project\UniProt_Query.py"

echo Running Expression IDs through UniProt
%pyPath% %scriptPath%

pause