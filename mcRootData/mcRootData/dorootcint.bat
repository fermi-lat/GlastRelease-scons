REM do rootcint for mc classes
cd ..\mcRootData
%ROOTSYS%\bin\rootcint -f mcRootDataCint.cxx -c RMcVertex.h RMcParticle.h RMcEvent.h mcRootDataLinkDef.h
