{
    new TBrowser;

// LOAD LIBRARIES:

    // Windows:
    gSystem->Load("c:/glast/software/digiRootData/v1r1/win32debug/digiRootData.dll");
    gSystem->Load("c:/glast/software/reconRootData/v1r2/win32debug/reconRootData.dll");
    gSystem->Load("c:/glast/software/mcRootData/v0r1/win32debug/mcRootData.dll");

    // Linux:  (for sun do:  /i386_linux/  ==>  /sun4x_56/  or  /sun4x_57/  etc.)
//    gSystem->Load("~/glast/software/digiRootData/v1r1/i386_linux/libdigiRootData.so");
//    gSystem->Load("~/glast/software/reconRootData/v1r2/i386_linux/libreconRootData.so");
//    gSystem->Load("~/glast/software/mcRootData/v0r1/i386_linux/libmcRootData.so");
}