{
    char* macros[]={
        "RootTreeAnalysis.cxx",
        "NtupleAnalysis.cxx",
        "chainTrees.cxx",
        "chainAll.cxx",
        "copyTree.cxx",
        "createEventList.cxx",
        "pruneTree.cxx"};
    char* libs[]={
        "mcRootData",
        "digiRootData",
        "reconRootData"};

    printf("\nLoading GLAST macros...\n");
    for(int i=0; i< sizeof(macros)/sizeof(char*); i++){
        printf("\t%s\n", macros[i]);
        gROOT->LoadMacro(macros[i]);
    }
    printf("\nLoading GLAST libraries...\n");
    // Load the appropriate version - depending on what OS 
    // we are running
    if (!strcmp(gSystem->GetName(), "WinNT")) {
        // nt version can use CMT for env vars.
        char buf[2000]; 
        for(int j=0; j< sizeof(libs)/sizeof(char*); j++){
            TString envStr(libs[j]);
            envStr.Append("ROOT");
            char *env = getenv(envStr.Data());
           // If we are not using CMT - we will use the ROOTANALYSIS env var to determine
           // where the GLAST ROOT libraries are located.
           char * raPath = getenv("ROOTANALYSIS");
           TString path;
           if (env) {
               path = env;
               path.Append("/");
               TString configStr(libs[j]);
               configStr.Append("CONFIG");
               char *config = getenv(configStr.Data());
               path.Append(config);
               path.Append("/");
               path.Append(libs[j]);
               path.Append(".dll");
           } else {
               path.Append(raPath);
               path.Append("/lib/");
               path.Append(libs[j]);
               path.Append(".dll");
           }
           gSystem->Load( path.Data() ); 
        }
       
    } else {  // UNIX
        gSystem->Load("libPhysics.so");
        gSystem->Load("libmcRootData.so");
        gSystem->Load("libdigiRootData.so");
        gSystem->Load("libreconRootData.so");
    }
    
}

