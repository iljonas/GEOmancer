@echo off
title UniProt Query Runner
set pyBase=C:\Users\%USERNAME%\AppData\Local\Programs\Python
for /F "tokens=*" %% in ('dir %pyBase% /A:D /O:D /T:C /B') do set pyVersion=%%i
set pyPath=%pyBase%\%pyVersion%\python.exe

if not exist %pyBase%\%pyVersion%\Lib\site-packages\Bio (
	%pyPath% -m pip install biopython)

C:\Users\Isaac\AppData\Local\Programs\Python\Python35\Lib\site-packages\Bio
pause