REM do rootcint for RECON ROOT classes, only if either Cint file
REM is missing
if not exist %RECONROOTDATAROOT%\reconRootData\reconRootDataCint.cxx goto :doit
if exist %RECONROOTDATAROOT%\reconRootData\reconRootDataCint.h goto :exit
:doit
REM Save the current directory, so we can return when finished
cd > dirstor.txt
copy %RECONROOTDATAROOT%\cmt\cd.cmd+dirstor.txt %RECONROOTDATAROOT%\changeDir.bat
del dirstor.txt

cd %RECONROOTDATAROOT%\reconRootData

%ROOT_PATH%\bin\rootcint -f reconRootDataCint.cxx -c -I%DIGIROOTDATAROOT% CalCluster.h CalRecon.h CalRangeRecData.h CalXtalRecData.h ReconEvent.h TkrCandHit.h TkrCandTrack.h TkrCovMat.h TkrFitHit.h TkrHitPlane.h TkrParams.h TkrRecon.h TkrSiCluster.h TkrTrack.h TkrVertex.h ReconLinkDef.h

REM Return to the original directory
call %RECONROOTDATAROOT%\changeDir.bat
del %RECONROOTDATAROOT%\changeDir.bat

:exit