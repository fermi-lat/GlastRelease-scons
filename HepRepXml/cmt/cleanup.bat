@echo off
if NOT DEFINED CMTROOT set CMTROOT=e:\Ric\glast\CMT\v1r14p20030620
call %CMTROOT%\mgr\setup.bat

set cmttempfile="%TEMP%\tmpsetup.bat"
%CMTROOT%\%CMTBIN%\cmt.exe -quiet cleanup -bat  -pack=HepRepXml -version=v0r3 -path=%~d0\Ric\glast\NewGlastRelease   %1 %2 %3 %4 %5 %6 %7 %8 %9 >%cmttempfile%
if exist %cmttempfile% call %cmttempfile%
if exist %cmttempfile% del %cmttempfile%
set cmttempfile=

