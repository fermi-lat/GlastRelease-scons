@echo off
if .%1==. (set tag=VC8 ) else set tag=%1
set CMTPATH=e:\Ric\glast\GlastRelease44;e:\Ric\glast\temp;e:\Ric\glast;
set tempfile=%HOME%\tmpsetup.bat
e:\Ric\glast/CMT/v1r16p20040701/VisualC/cmt.exe -quiet -bat -pack=HepRepXml -version=v0r7 setup -tag=%tag% >"%tempfile%"
if exist "%tempfile%" call "%tempfile%"
if exist "%tempfile%" del "%tempfile%"
set PATH=%LD_LIBRARY_PATH%;%PATH%