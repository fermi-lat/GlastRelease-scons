REM $Header$
@echo off

setlocal
set PYTHONPATH=%CALIBGENCALROOT%\python\lib;%ROOTSYS%\bin;%PYTHONPATH%;
python %CALIBGENCALROOT%\python\checkXML.py %1 %2 %3 %4 %5 %6 %7 %8 %9
endlocal

