set SYSTEST=ACDDigi
call ..\..\..\cmt\setup.bat
call %GLEAMROOT%\src\systests\%SYSTEST%\noiseTest.bat
del *-mc.root *-digi.root *-recon.root
