@echo off
if .%1==. (set tag=Win32Debug ) else set tag=%1
set tempfile=%HOMEDRIVE%%HOMEPATH%tmpsetup.bat
c:\glast\CMT\v1r10\VisualC\cmt.exe -quiet -bat -pack=G4Generator -version=v0r3 setup -tag=%tag% >%tempfile%
if exist %tempfile% call %tempfile%
if exist %tempfile% del %tempfile%
set PATH=%LD_LIBRARY_PATH%;%PATH%