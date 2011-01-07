
#include <reconRootData/ReconEvent.h>
#include <commonRootData/RootObj.h>
#include <commonRootData/RootDataUtil.h>
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom.h"
#include "TCollection.h"
#include "TVector3.h"
#include "Riostream.h"
#include <string>
#include "facilities/Util.h"


/** @file testReconClasses.cxx
* @brief This defines a test routine for the Reconstruction ROOT classes.
*
* This program creates a new Recon ROOT file, and the opens it up again
* for reading.  The contents are checked against the values known to be stored
* during the original writing of the file.
* The contents of the file are printed to the screen.
* The program returns 0 if the test passed.
* If failure, the program returns -1.
*
* $Header$
*/

const UInt_t RUN_NUM = 1;

int checkReconEvent(ReconEvent *evt, UInt_t ievent) {
    
    if (evt->getRunId() != RUN_NUM) {
        std::cout << "Run Id is wrong: " << evt->getRunId() << std::endl;
        return -1;
    }
    if (evt->getEventId() != ievent) {
        std::cout << "Event Id is wrong: " << evt->getEventId() << std::endl;
        return -1;
    }

    return 0;
    
}

/// Read in the ROOT file just generated via the write method
int read(const char* fileName, int numEvents) {
    
    TFile *f = new TFile(fileName, "READ");

    AcdRecon::fixAcdStreamer(f->GetVersion()) ;

   
    TObjArray * randomNumbers = (TObjArray*)f->Get("randomNumbers") ;
   
    TTree *t = (TTree*)f->Get("Recon");
    ReconEvent *evt = 0;
    t->SetBranchAddress("ReconEvent", &evt);

    std::cout << "Opened the ROOT file for reading" << std::endl;

    UInt_t ievent;
    for (ievent = 0; ievent < (UInt_t) numEvents; ievent++) {
        
        Double_t randomNumber = ((RootObj<Double_t> *)randomNumbers->At(ievent))->value() ;
        
        t->GetEvent(ievent);
        std::cout << "ReconEvent ievent = " << ievent << std::endl;
        evt->Print() ;
        
        if (checkReconEvent(evt, ievent) < 0)
         { return -1 ; }
        
        // prepare data
        AcdRecon * acd = evt->getAcdRecon() ;
        CalRecon * cal = evt->getCalRecon() ;
        TkrRecon * tkr = evt->getTkrRecon() ;
        reconRootData::AdfRecon * adf = evt->getAdfRecon();
        if ((!acd)||(!cal)||(!tkr)||(!adf))
         { return -1 ; }
        acd->Print() ;
        cal->Print();
        tkr->Print();
        adf->Print();
        AcdRecon acdRef ;
        CalRecon calRef ;
        TkrRecon tkrRef ;
        reconRootData::AdfRecon adfRef;
         
        // usual tests
        acdRef.Fake(ievent,randomNumber) ;
        if (!acd->CompareInRange(acdRef))
         { return -1 ; }
        calRef.Fake(ievent,randomNumber) ;
        if (!cal->CompareInRange(calRef))
         { return -1 ; }
        tkrRef.Fake(ievent,randomNumber) ;
        if (!tkr->CompareInRange(tkrRef))
         { return -1 ; }
        adfRef.Fake(ievent,randomNumber);
        if (!adf->CompareInRange(adfRef))
         { return -1; }
         
        // opposite tests
        std::cout<<"===== on purpose errors ====="<<std::endl ;
        acdRef.Clear() ;
        acdRef.Fake(ievent,randomNumber*2.) ;
        if (acd->CompareInRange(acdRef))
         { return -1 ; }
        calRef.Clear() ;
        calRef.Fake(ievent,randomNumber*2.) ;
        if (cal->CompareInRange(calRef))
         { return -1 ; }
        tkrRef.Clear() ;
        tkrRef.Fake(ievent,randomNumber*2.) ;
        if (tkr->CompareInRange(tkrRef))
         { return -1 ; }
        adfRef.Clear();
        adfRef.Fake(ievent,randomNumber*2.);
        if (!adf->CompareInRange(adfRef))
         { return -1; }
        std::cout<<"============================="<<std::endl ;
        
        evt->Clear();
    }

    f->Close();
    delete f;

    return 0;
}

/// Create a new ROOT file
int write(const char* fileName, int numEvents) {
    
    Int_t buffer = 64000 ;
    Int_t splitLevel = 1 ;

    TFile * f =  new TFile(fileName,"RECREATE") ;
    
    TTree * t = new TTree("Recon","Recon") ;
    ReconEvent * ev = new ReconEvent() ;
    t->Branch("ReconEvent", "ReconEvent", &ev, buffer, splitLevel);
    std::cout << "Created new ROOT file" << std::endl;

    TRandom randGen ;
    TObjArray randomNumbers ;
    
    Int_t ievent ;
    for (ievent = 0; ievent < numEvents; ievent++) {

        RootObj<Double_t> * randomNumber = new RootObj<Double_t>(randGen.Rndm()) ;
        randomNumbers.Add(randomNumber) ;

        AcdRecon * acdRec = new AcdRecon() ;
        acdRec->Fake(ievent,randomNumber->value()) ;

        AcdReconV2 * acdRecV2 = new AcdReconV2() ;
        acdRecV2->Fake(ievent,randomNumber->value()) ;

        CalRecon * calRec = new CalRecon() ;
        calRec->Fake(ievent,randomNumber->value()) ;

        TkrRecon * tkrRec = new TkrRecon();
        tkrRec->Fake(ievent,randomNumber->value()) ;        

        reconRootData::AdfRecon * adfRec = new reconRootData::AdfRecon();
        adfRec->Fake(ievent,randomNumber->value());

        ev->initialize(ievent, RUN_NUM, tkrRec, calRec, acdRec, acdRecV2);
        ev->initAdf(adfRec);
        ev->initEventFlags(1) ;
        
        t->Fill();
        ev->Clear();
    }

    std::cout << "Filled ROOT file with " << numEvents << " events" << std::endl;
    delete ev ;
    
    f->cd() ;
    randomNumbers.Write("randomNumbers",TObject::kSingleKey) ;
    
    f->Write() ;
    f->Close();
    delete f;
    
    randomNumbers.Delete() ;
    randomNumbers.Clear() ;
    
    return(0);
}


/// Main program
/// Return 0 for success.
/// Returns -1 for failure.
int main(int argc, char **argv) {
    
    std::string fileStr("$(TESTDATAROOT)/recon-v9r25-v5r26.root");
    facilities::Util::expandEnvVar(&fileStr);
    //char * fileName = "recon-v9r11-v5r10.root" ;
    const char *fileName = fileStr.c_str();
    int numEvents = 10 ;
    
    TString writeString("write") ;
    bool writeMode = false ;

    int n = 1;
    if (argc > 1) {
      if (writeString.CompareTo(argv[n++])==0)
        { writeMode = true ; }
    } 
    if (argc > 2) {
        fileName = argv[n++];
    } 
    if (argc > 3) {
        numEvents = atoi(argv[n++]);
    } 

    int sc = 0 ;
    try 
     {
      if (writeMode)
       { sc = write(fileName,numEvents) ; }
      else
       { sc = read(fileName,numEvents) ; }
     }
    catch (...)
     {
      std::cout<<"AN UNKNOWN EXCEPTION HAS BEEN RAISED"<<std::endl ;
      sc = 1 ;
     }
     
     
    if (sc == 0)
     {
      if (writeMode)
       { std::cout << "RECON ROOT file writing succeeded!" << std::endl ; }
      else
       { std::cout << "RECON ROOT file reading succeeded!" << std::endl ; }
     }
    else
     {
      if (writeMode)
       { std::cout << "FAILED recon writing" << std::endl ; }
      else
       { std::cout << "FAILED recon reading" << std::endl ; }
     }

    return(sc);
}
