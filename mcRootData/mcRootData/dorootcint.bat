REM do rootcint for mc classes
cd ..\mcRootData
%ROOTSYS%\bin\rootcint -f mcRootDataCint.cxx -c McVertex.h McParticle.h McEvent.h mcRootDataLinkDef.h
