@echo off
if .%1==. (set tag=Win32Debug ) else set tag=%1
set tempfile=%HOME%\tmpsetup.bat
c:\Programs\CMT\v1r12p20021129\VisualC\cmt.exe -quiet -bat -pack=AnalysisNtuple -version=v0 setup -tag=%tag% >"%tempfile%"
if exist "%tempfile%" call "%tempfile%"
if exist "%tempfile%" del "%tempfile%"
set PATH=%LD_LIBRARY_PATH%;%PATH%