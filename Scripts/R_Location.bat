@ECHO off
SET Rbase="C:\Program Files\R"
FOR /F "tokens=*" %%i IN ('DIR %Rbase% /A:D /O:D /T:C /B') DO SET Rversion=%%i
SET Rpath="%Rbase:"=%\%Rversion%\bin\Rscript.exe"
ECHO %Rpath%