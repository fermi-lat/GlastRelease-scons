REM do rootcint for MC ROOT classes, only if either Cint file
REM is missing
if not exist %MCROOTDATAROOT%\mcRootData\mcRootDataCint.cxx goto :doit
if exist %MCROOTDATAROOT%\mcRootData\mcRootDataCint.h goto :exit
:doit
REM Save the current directory, so we can return when finished
cd > dirstor.txt
copy %DIGIROOTDATAROOT%\cmt\cd.cmd+dirstor.txt %DIGIROOTDATAROOT%\changeDir.bat
del dirstor.txt

cd %DIGIROOTDATAROOT%\digiRootData
REM Remove old cint files
del EventCint.cxx
del EventCint.h

%ROOT_PATH%\bin\rootcint -f digiRootDataCint.cxx -c CalDigi.h CalXtalId.h CalXtalReadout.h DigiEvent.h EventLinkDef.h

REM Return to the original directory
call %DIGIROOTDATAROOT%\changeDir.bat
del %DIGIROOTDATAROOT%\changeDir.bat
