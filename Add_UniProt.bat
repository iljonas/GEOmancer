@echo off
TITLE Add_UniProt

C:\Users\%USERNAME%\Documents\GEOmancer\Scripts\R_Location.bat > tempFile
SET /P Rpath= < tempFile

C:\Users\%USERNAME%\Documents\GEOmancer\Scripts\Python_Location.bat > tempFile
SET /P pyPath= < tempFile

DEL tempFile

SET sourcePath="C:\Users\%USERNAME%\Documents\GEOmancer\Scripts
SET queryPath=%sourcePath%\UniProt_Query.py"
SET transformPath=%sourcePath%\UniProt_Transformer.R"

SET /P input="Enter source file name (do not include file extension): "

ECHO Running Expression IDs through UniProt
%pyPath% %queryPath% %input%

ECHO.
ECHO Combining UniProt output with project master file
%RPath% %transformPath% %input%

PAUSE