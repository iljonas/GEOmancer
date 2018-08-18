@ECHO off
TITLE Chart_Builder

C:\Users\%USERNAME%\Documents\GEOmancer\Scripts\Python_Location.bat > tempFile
PAUSE
SET /P pyPath= < tempFile
DEL tempFile

SET /P master="Enter master file name (do not include file extension): "
SET scriptPath=C:\Users\%USERNAME%\Documents\GEOmancer\Scripts\Auto_Query.py

%pyPath% %scriptPath%

PAUSE