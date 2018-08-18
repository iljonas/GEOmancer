@ECHO off
SET pyBase=C:\Users\%USERNAME%\AppData\Local\Programs\Python
FOR /F "tokens=*" %%i IN ('DIR %pyBase% /A:D /O:D /T:C /B') DO SET pyVersion=%%i
SET pyPath=%pyBase%\%pyVersion%\python.exe
ECHO %pyPath%