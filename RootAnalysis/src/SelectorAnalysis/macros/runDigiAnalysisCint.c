{
    // Run setup.bat before running
    #include "Riostream.h"
    gROOT->Reset();    // roll back CINT context to last Save

    gSystem->Load("libmcRootData.so");
    gSystem->Load("libdigiRootData.so");
    gSystem->Load("libreconRootData.so");
    gSystem->Load("libglastSelector.so");

    // Create the GLAST folders
    TFolder *glastFolder = gROOT->GetRootFolder()->AddFolder("glast","glast top level folders");
    TFolder *histoFolder = glastFolder->AddFolder("histograms", "histogram folder");
    TFolder *taskFolder = glastFolder->AddFolder("tasks", "task folder");


    // Setup our tasks
    GlastDigiTask *digiTask = new GlastDigiTask("DigiTask", "Process digis");
    AcdDigiTask *acdT = new AcdDigiTask("AcdDigiTask", "Process one set of acd digis from one event");
    // Make AcdTask a subtask of the main GlastDigiTask
    digiTask->Add(acdT);
    taskFolder->Add(digiTask);

    // Now open our data file and run the selector
    TFile f("digi_1500.root", "READ");
    TTree *t = (TTree*)f.Get("Digi");
    // Process 50 events starting with event 0
    DigiSelector *selector = new DigiSelector();
    t->Process(selector,"",500,0);    
    TBrowser tb;

}







