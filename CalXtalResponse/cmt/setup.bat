rem Setting CalXtalResponse v0r0 in %~d0%~p0
@echo off
if NOT DEFINED CMTROOT set CMTROOT=d:\GLAST\CMT\CMT\v1r14p20031120 & set PATH=%CMTROOT%\%CMTBIN%;%PATH% & set CMTBIN=VisualC & if not defined CMTCONFIG set CMTCONFIG=%CMTBIN%

set cmttempfile="%TEMP%\tmpsetup.bat"
%CMTROOT%\%CMTBIN%\cmt.exe -quiet setup -bat  -pack=CalXtalResponse -version=v0r0 -path=%~d0%~p0..\..\..   %1 %2 %3 %4 %5 %6 %7 %8 %9 >%cmttempfile%
if exist %cmttempfile% call %cmttempfile%
if exist %cmttempfile% del %cmttempfile%
set cmttempfile=

