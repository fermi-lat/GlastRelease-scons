@echo off
REM $Header$
REM do rootcint for MC ROOT classes, only if both Cint classes exist
if not exist %MCROOTDATAROOT%\mcRootData\mcRootDataCint.cxx goto :doit
if exist %MCROOTDATAROOT%\mcRootData\mcRootDataCint.h goto :exit
:doit
REM Save the current directory, so we can return when finished
cd > dirstor.txt
copy %MCROOTDATAROOT%\cmt\cd.cmd+dirstor.txt %MCROOTDATAROOT%\changeDir.bat
del dirstor.txt
cd %MCROOTDATAROOT%\mcRootData

REM Remove old cint files
if  exist  mcRootDataCint.cxx del mcRootDataCint.cxx
if  exist  mcRootDataCint.h del mcRootDataCint.h

%ROOT_PATH%\bin\rootcint -f mcRootDataCint.cxx -c McEvent.h McIntegratingHit.h McParticle.h McPositionHit.h VolumeIdentifier.h mcRootDataLinkDef.h 

REM Return to the original directory
call %MCROOTDATAROOT%\changeDir.bat
del %MCROOTDATAROOT%\changeDir.bat
:exit