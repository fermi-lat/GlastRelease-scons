REM $Header$
@echo off

if not defined CALIBGENCALROOT goto :ERROR

setlocal
set PYTHONPATH=%CALIBGENCALROOT%\python\lib;%PYTHONPATH%;
python %CALIBGENCALROOT%\python\muSlopeXML2TXT.py %1 %2 %3 %4 %5 %6 %7 %8 %9
endlocal
goto :EXIT

:ERROR
echo ERROR: muSlopeXML2TXT: CALIBGENCALROOT must be defined

:EXIT
