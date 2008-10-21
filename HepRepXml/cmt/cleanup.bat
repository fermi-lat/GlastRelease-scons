@echo off
if NOT DEFINED CMTROOT set CMTROOT=c:\GlastPrograms\CMT\v1r18p20061003& set PATH=%CMTROOT%\%CMTBIN%;%PATH%& set CMTBIN=VisualC& if not defined CMTCONFIG set CMTCONFIG=%CMTBIN%

set cmttempfile="%TEMP%\tmpsetup.bat"
%CMTROOT%\%CMTBIN%\cmt.exe cleanup -bat  -pack=HepRepXml -version=v0r8p2 -path="%~d0%~p0..\..\.."   %1 %2 %3 %4 %5 %6 %7 %8 %9 >%cmttempfile%
if exist %cmttempfile% call %cmttempfile%
if exist %cmttempfile% del %cmttempfile%
set cmttempfile=

