REM do rootcint for Digi ROOT classes

REM Save the current directory, so we can return when finished
cd > dirstor.txt
copy %DIGIROOTDATAROOT%\cmt\cd.cmd+dirstor.txt %DIGIROOTDATAROOT%\changeDir.bat
del dirstor.txt

cd %DIGIROOTDATAROOT%\digiRootData
REM Remove old cint files
del EventCint.cxx
del EventCint.h

%ROOT_PATH%\bin\rootcint -f EventCint.cxx -c AcdHeader.h AcdId.h AcdTile.h CalHeader.h CalDigi.h CalLogId.h CalLogReadout.h CalLog.h ESAPID.h DigiEvent.h L1T.h LiveTime.h LogId.h StripId.h SubSystemHeader.h Tagger.h TkrHeader.h TkrLayer.h EventLinkDef.h

REM Return to the original directory
call %DIGIROOTDATAROOT%\changeDir.bat
del %DIGIROOTDATAROOT%\changeDir.bat
