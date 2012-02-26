#include "Riostream.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom.h"
#include "TCollection.h"
#include "gcrSelectRootData/GcrSelectEvent.h"
#include <cstdlib>

/** @file testGcrClasses.cxx
* @brief This defines a test routine for the GCR ROOT classes.
*
* This program create a new Gcr ROOT file, and the opens it up again
* for reading.  The contents are checked against the values known to be stored
* during the original writing of the file.
* The contents of the file are printed to the screen.
* The program returns 0 if the test passed.
* If failure, the program returns -1.
*
* $Header$
*/
const UInt_t runNum = 1;
Float_t randNum;


int checkGcrSelectEvent(GcrSelectEvent *evt, UInt_t ievent) {

    GcrSelectEvent evtRef;
    evtRef.Fake(ievent, randNum);
    return evt->CompareInRange(evtRef, "GcrSelectEvent");

}



/// Read in the ROOT file just generated via the write method
int read(const char* fileName, int numEvents) {
    TFile *f = new TFile(fileName, "READ");
    TTree *t = (TTree*)f->Get("GcrSelect");
    GcrSelectEvent *evt = 0;
    t->SetBranchAddress("GcrSelectEvent", &evt);
    
    std::cout << "Opened the ROOT file for reading" << std::endl;
    
    Int_t ievent;
    for (ievent = 0; ievent < numEvents; ievent++) {
        t->GetEvent(ievent);
        std::cout << "GcrSelectEvent ievent = " << ievent << std::endl;
        evt->Print();
        if (checkGcrSelectEvent(evt, ievent) < 0) return -1;

    }
    
    f->Close();
    delete f;
    
    return 0;
}

/// Create a new Monte Carlo ROOT file
int write(const char* fileName, int numEvents) {
    Int_t buffer = 64000;
    Int_t splitLevel = 1;
    
    TFile *f =  new TFile(fileName, "RECREATE");
    TTree *t = new TTree("GcrSelect", "GcrSelect");
    GcrSelectEvent *ev = new GcrSelectEvent();
    t->Branch("GcrSelectEvent", "GcrSelectEvent", &ev, buffer, splitLevel);
    
    std::cout << "Created new ROOT file" << std::endl;
    
    TRandom randGen;
    Int_t ievent;
    randNum = randGen.Rndm();

    for (ievent = 0; ievent < numEvents; ievent++) {
        
        ev->Fake(ievent, randNum);
        ev->Print();

        t->Fill();
        ev->Clear();
    }
    
    std::cout << "Filled ROOT file with " << numEvents << " events" << std::endl;
    delete ev;
    
    f->Write();
    f->Close();
    delete f;
    return(0);
}


/// Main program
/// Return 0 for success.
/// Returns -1 for failure.
int main(int argc, char **argv) {
    const char *fileName = "gcr.root";
    int n = 1;
    int numEvents = 10;
    if (argc > 1) {
        fileName = argv[n++];
    } 
    if (argc > 2) {
        numEvents = atoi(argv[n++]);
    } 
    
    int sc = 0;
    sc = write(fileName, numEvents);
    sc = read(fileName, numEvents);
    
    if (sc == 0) {
        std::cout << "GCR ROOT file writing and reading succeeded!" << std::endl;
    } else {
        std::cout << "FAILED GCR writing and reading" << std::endl;
    }
    
    return(sc);
}


