@ECHO off
TITLE Check for BioPython

C:\Users\%USERNAME%\Documents\GEOmancer\Scripts\Python_Location.bat > tempFile
SET /P pyPath= < tempFile
DEL tempFile

:START
IF NOT EXIST %pyBase%\%pyVersion%\Lib\site-packages\Bio (
	SET /P downloadResponse="This BLAST query method requires BioPython. Do you like to install it (y/n)? "
)
IF /I %downloadResponse% NEQ y IF /I %downloadResponse% NEQ n GOTO START

IF /I %downloadResponse% EQU y (
	%pyPath% -m pip install biopython %*
)

PAUSE