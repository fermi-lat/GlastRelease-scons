REM do rootcint for DIGI ROOT classes, only if either Cint file
REM is missing
if not exist %DIGIROOTDATAROOT%\digiRootData\digiRootDataCint.cxx goto :doit
if exist %DIGIROOTDATAROOT%\digiRootData\digiRootDataCint.h goto :exit
:doit
REM Save the current directory, so we can return when finished
cd > dirstor.txt
copy %DIGIROOTDATAROOT%\cmt\cd.cmd+dirstor.txt %DIGIROOTDATAROOT%\changeDir.bat
del dirstor.txt

cd %DIGIROOTDATAROOT%\digiRootData

%ROOT_PATH%\bin\rootcint -f digiRootDataCint.cxx -c AcdDigi.h AcdId.h CalDigi.h CalXtalId.h CalXtalReadout.h DigiEvent.h GlastAxis.h TkrDigi.h TowerId.h EventLinkDef.h

REM Return to the original directory
call %DIGIROOTDATAROOT%\changeDir.bat
del %DIGIROOTDATAROOT%\changeDir.bat
:exit
