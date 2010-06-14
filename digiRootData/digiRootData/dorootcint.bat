REM do rootcint for Event classes
cd ..\digiRootData
%ROOTSYS%\bin\rootcint -f EventCint.cxx -c ACDTile.h CalHit.h ESAPID.h Event.h L1T.h LogID.h StripID.h Tagger.h TileID.h TkrLayer.h EventLinkDef.h
