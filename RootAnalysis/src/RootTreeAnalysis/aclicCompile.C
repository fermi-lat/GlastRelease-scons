{
// source setup in cmt directory first, then start ROOT and run this macro
gROOT->Reset();    // roll back CINT context to last Save
gSystem->Load("$COMMONROOTDATAROOT/$CMTCONFIG/libcommonRootData.so");
gSystem->Load("$MCROOTDATAROOT/$CMTCONFIG/libmcRootData.so");
gSystem->Load("$DIGIROOTDATAROOT/$CMTCONFIG/libdigiRootData.so");
gSystem->Load("$RECONROOTDATAROOT/$CMTCONFIG/libreconRootData.so");
gSystem->AddIncludePath(" -I$ENUMSROOT" );
gSystem->AddIncludePath(" -I$COMMONROOTDATAROOT" );
gSystem->AddIncludePath(" -I$MCROOTDATAROOT" );
gSystem->AddIncludePath(" -I$DIGIROOTDATAROOT" );
gSystem->AddIncludePath(" -I$RECONROOTDATAROOT" );
gROOT->LoadMacro("RootTreeAnalysis.cxx++");
}
