REM do rootcint for Event classes
cd ..\digiRootData
%ROOTSYS%\bin\rootcint -f EventCint.cxx -c AcdHeader.h AcdId.h AcdTile.h CalHeader.h CalLog.h CalLogId.h ESAPID.h Event.h L1T.h StripId.h SubSystemHeader.h Tagger.h TkrHeader.h TkrLayer.h EventLinkDef.h
