@ECHO off
TITLE Chart_Builder

C:\Users\%USERNAME%\Documents\GEOmancer\Scripts\R_Location.bat > tempFile
SET /P Rpath= < tempFile

SET scriptPath="C:\Users\%USERNAME%\Documents\GEOmancer\Scripts\GEO_Assembler.R"
SET /P series="Enter Series ID: "
SET /P master="Enter master file name (do not include file extension): "
SET /P exclusion="Enter column exclusions (no extra spaces): "

ECHO Analyzing Series and Platform Files...
%Rpath% %scriptPath% %series% %master% %exclusion%

PAUSE