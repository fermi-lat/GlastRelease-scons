if NOT DEFINED CMTROOT set CMTROOT=e:\Ric\glast\CMT\v1r14
call %CMTROOT%\mgr\setup.bat
set tempfile="%HOMEDRIVE%%HOMEPATH%\tmpsetup.bat"
%CMTROOT%\%CMTBIN%\cmt.exe -quiet cleanup -bat -path=e:\Ric\glast\GlastRelease %1 %2 %3 %4 %5 %6 %7 %8 %9 >%tempfile%
if exist %tempfile% call %tempfile%
if exist %tempfile% del %tempfile%

