REM do rootcint for Event classes
cd ..\digiRootData
%ROOTSYS%\bin\rootcint -f EventCint.cxx -c AcdHeader.h AcdId.h AcdDigi.h CalHeader.h CalDigi.h CalLogId.h ESAPID.h Event.h L1T.h StripId.h SubSystemHeader.h Tagger.h TkrHeader.h TkrDigi.h EventLinkDef.h
