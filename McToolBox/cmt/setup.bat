@echo off
if .%1==. (set tag=VC8debug ) else set tag=%1
set CMTPATH=C:\Glast\Development;C:\Glast\GlastRelease_v4r4;C:\Glast\RootAnalysis;C:\Glast\tools\Glast_external\Gaudi\v12r0p0;C:\Glast\tools;
set tempfile=%HOME%\tmpsetup.bat
C:\Glast\tools/CMT/v1r16p20040701/VisualC/cmt.exe -quiet -bat -pack=McToolBox -version=v1r0p0 setup -tag=%tag% >"%tempfile%"
if exist "%tempfile%" call "%tempfile%"
if exist "%tempfile%" del "%tempfile%"
set PATH=%LD_LIBRARY_PATH%;%PATH%