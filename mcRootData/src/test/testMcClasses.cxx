
#include <mcRootData/McEvent.h>
#include <commonRootData/RootDataUtil.h>
#include "Riostream.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom.h"

/** @file testMcClasses.cxx
* @brief This defines a test routine for the Monte Carlo ROOT classes.
*
* This program create a new Monte Carlo ROOT file, and the opens it up again
* for reading.  The contents are checked against the values known to be stored
* during the original writing of the file.
* The contents of the file are printed to the screen.
* The program returns 0 if the test passed.
* If failure, the program returns -1.
*
* $Header$
*/

const UInt_t RUN_NUM = 1 ;
Float_t RAND_NUM ;

int read( char * fileName, unsigned int numEvents) {
    // Purpose and Method: Read in the ROOT file just generated via the 
    //   write method

    TFile *f = new TFile(fileName, "READ");
    TTree *t = (TTree*)f->Get("Mc");
    McEvent *evt = 0;
    t->SetBranchAddress("McEvent", &evt);
    
    std::cout << "Opened the ROOT file for reading" << std::endl;
    
    UInt_t iEvent ;
    for ( iEvent = 0 ; iEvent < numEvents ; ++iEvent ) {
        
        t->GetEvent(iEvent);
     
        std::cout << "McEvent iEvent = " << iEvent << std::endl;
        evt->Print() ;
                
        // DC: I CANNOT MAKE A McEvent::Fake() and
        // a McEvent::CompareInRange(), because McEvent
        // is designed as a singleton.
        
        if (!evt->CompareToFake(iEvent,RUN_NUM,RAND_NUM)) {
          return -1 ;
        }
        
    }
    
    f->Close();
    delete f;
    
    return(0);
}

/// Create a new ROOT file
int write( char * fileName, UInt_t numEvents ) {

    Int_t buffer = 64000;
    Int_t splitLevel = 1;
    
    TFile *f =  new TFile(fileName, "RECREATE");
    TTree *t = new TTree("Mc", "Mc");
    McEvent * evt = new McEvent() ;
    t->Branch("McEvent", "McEvent", &evt, buffer, splitLevel);
    std::cout << "Created new ROOT file" << std::endl;
    
    UInt_t iEvent ;
    for (iEvent = 0; iEvent < numEvents; iEvent++) {
        evt->Fake(iEvent,RUN_NUM,RAND_NUM) ;
        t->Fill() ;
    }
    
    std::cout << "Filled ROOT file with " << numEvents << " events" << std::endl;
    delete evt ;
    
    f->Write();
    f->Close();
    delete f;
    return(0);
}


/// Main program
/// Return 0 for success.
/// Returns -1 for failure.
int main(int argc, char **argv) {

    char *fileName = "mc.root";
    int n =1 ;
    unsigned int numEvents =10 ;
    if (argc > 1) {
        fileName = argv[n++];
    } 
    if (argc > 2) {
        numEvents = atoi(argv[n++]);
    } 
    
    TRandom randGen ;
    RAND_NUM = randGen.Rndm() ;

    int sc = 0;
    try 
     {
      sc = write(fileName, numEvents);
      sc = read(fileName, numEvents);
     }
    catch (...)
     {
      std::cout<<"AN UNKNOWN EXCEPTION HAS BEEN RAISED"<<std::endl ;
      sc = 1 ;
     }
    
    if (sc == 0) {
        std::cout << "MC ROOT file writing and reading suceeded!" << std::endl;
    } else {
        std::cout << "FAILED" << std::endl;
    }
    
    return(sc);
}


