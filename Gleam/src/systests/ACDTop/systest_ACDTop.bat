set SYSTEST=ACDTop
call ..\..\..\cmt\setup.bat
call %GLEAMROOT%\src\systests\%SYSTEST%\test_Mip.bat

# make sure Gleam uses local job options file
..\..\..\Win32Debug\Gleam.exe %GLEAMROOT%\src\systests\%SYSTEST%\muon_300MeV_p
atch_run1.txt
..\..\..\Win32Debug\Gleam.exe %GLEAMROOT%\src\systests\%SYSTEST%\muon_300MeV_p
atch_run2.txt
..\..\..\Win32Debug\Gleam.exe %GLEAMROOT%\src\systests\%SYSTEST%\muon_300MeV_p
atch_run3.txt
..\..\..\Win32Debug\Gleam.exe %GLEAMROOT%\src\systests\%SYSTEST%\muon_300MeV_p
atch_run4.txt
..\..\..\Win32Debug\Gleam.exe %GLEAMROOT%\src\systests\%SYSTEST%\muon_300MeV_p
atch_run5.txt

# RUN root macro contains boilerplate calls for
# generating root histograms from Gleam output
root -b -q RUN_Windows

mv Histograms.root ${SYSTEST}Histos.root

rm *mc.root *digi.root *recon.root

