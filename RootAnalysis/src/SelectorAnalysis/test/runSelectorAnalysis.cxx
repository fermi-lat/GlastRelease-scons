#include "TROOT.h"
#include "TSystem.h"
#include "TBrowser.h"
#include <iostream>
#include <string>
#include "GlastDigiTask.h"
#include "AcdDigiTask.h"
#include "DigiSelector.h"

/** @file runTreeAnalysis
@brief Simple test program to exercise RootTreeAnalysis
*/
int main(int argn, char** argc) {
    
#ifdef WIN32
  //  gSystem->Load("libTree.dll");
#endif
    unsigned int numEvents = 100;
    const char* path = ::getenv("ROOTANALYSISROOT");
    std::string digiFileName(path);
    digiFileName += "/src/test/digi.root";
    if (argn > 1) digiFileName = argc[1];
    std::cout << "Setup for Processing" << std::endl;

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
    TFile f(digiFileName.c_str(), "READ");
    TTree *t = (TTree*)f.Get("Digi");
    // Process 50 events starting with event 0
    DigiSelector *selector = new DigiSelector();
    t->Process(selector,"",50,0);    
    std::cout << "Done Processing " << std::endl;
    TBrowser *tb = new TBrowser("Browser","NewBrowser");

    return 0;
}







