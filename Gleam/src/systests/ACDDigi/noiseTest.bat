set SYSTEST=ACDDigi
call ..\..\..\cmt\setup.bat
%GLEAMROOT%\Win32Debug\Gleam.exe %GLEAMROOT%\src\systests\%SYSTEST%\gamma_100MeV_miss_0stddev.txt
%GLEAMROOT%\Win32Debug\Gleam.exe %GLEAMROOT%\src\systests\%SYSTEST%\gamma_100MeV_miss_0.01stddev.txt
%GLEAMROOT%\Win32Debug\Gleam.exe %GLEAMROOT%\src\systests\%SYSTEST%\gamma_100MeV_miss_0.02stddev.txt
%GLEAMROOT%\Win32Debug\Gleam.exe %GLEAMROOT%\src\systests\%SYSTEST%\gamma_100MeV_miss_0.04stddev.txt
%GLEAMROOT%\Win32Debug\Gleam.exe %GLEAMROOT%\src\systests\%SYSTEST%\gamma_100MeV_miss_0.1stddev.txt
%ROOTSYS%\bin\root -b -q RUN_noise_Windows

