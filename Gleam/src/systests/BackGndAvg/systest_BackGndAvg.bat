set SYSTEST=BackGndAvg
call ..\..\..\cmt\setup.bat
%GLEAMROOT%\Win32Debug\Gleam.exe %GLEAMROOT%\src\systests\%SYSTEST%\jobOptions.txt
%ROOTSYS%\bin\root -b -q RUN_Windows
move Histograms.root %SYSTEST%Histos.root
del mc.root digi.root recon.root