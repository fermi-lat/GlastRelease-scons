REM do rootcint for mc classes
cd ..\mcRootData
%ROOTSYS%\bin\rootcint -f mcRootDataCint.cxx -c MCVertex.h MCParticle.h mcRootDataLinkDef.h
