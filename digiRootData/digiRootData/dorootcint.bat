REM do rootcint for Event classes
cd ..\digiRootData
%ROOTSYS%\bin\rootcint -f EventCint.cxx -c AcdHeader.h AcdTile.h CalHeader.h CalHit.h ESAPID.h Event.h L1T.h LogID.h StripID.h Tagger.h TileID.h TkrHeader.h TkrLayer.h EventLinkDef.h
