REM do rootcint for Recon classes
cd ..\reconRootData
%ROOTSYS%\bin\rootcint -f ReconCint.cxx -c -I%DIGIROOTDATAROOT% CalCluster.h CalLogEne.h CalRecon.h Recon.h ReconHeader.h TkrGamma.h TkrHit.h TkrLocator.h TkrRecon.h TkrSiCluster.h TkrTrack.h ReconLinkDef.h
