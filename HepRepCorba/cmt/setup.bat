@echo off
if .%1==. (set tag=VC7debug ) else set tag=%1
set tempfile=%HOME%\tmpsetup.bat
e:\Ric\glast\CMT\v1r14\VisualC\cmt.exe -quiet -bat -pack=HepRepCorba -version=v0r1 setup -tag=%tag% >"%tempfile%"
if exist "%tempfile%" call "%tempfile%"
if exist "%tempfile%" del "%tempfile%"
set PATH=%LD_LIBRARY_PATH%;%PATH%