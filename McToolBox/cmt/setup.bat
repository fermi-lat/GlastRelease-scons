@echo off
if .%1==. (set tag=VC8debug ) else set tag=%1
set tempfile=%HOME%\tmpsetup.bat
C:\GLAST\tools\CMT\v1r14p20031120\VisualC\cmt.exe -quiet -bat -pack=McToolBox -version=v0r0 setup -tag=%tag% >"%tempfile%"
if exist "%tempfile%" call "%tempfile%"
if exist "%tempfile%" del "%tempfile%"
set PATH=%LD_LIBRARY_PATH%;%PATH%