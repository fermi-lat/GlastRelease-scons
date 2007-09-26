/** 
*@class afterglow
*
*@brief Root Tree Ananlysis 
*
* This class is for Trigger Root Tree Ananlysis of the afterglow in the TKR
*
* Version 0.1 24-August-2006 Martin Kocian Creation
*/

#ifndef AFTERGLOW_H
#define AFTERGLOW_H 1

// Need these includes if we wish to compile this code
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TMath.h>
#include <TStyle.h>
#include <TH2D.h>
#include <TProfile.h>
#include "../calibGenTRG/RootTreeAnalysis.h"
#include <TH1D.h>

#include  <iostream> 
#include  <iomanip> 
#include  <fstream> 

class afterglow : public RootTreeAnalysis {

public :

  afterglow(); 

  virtual ~afterglow() {;}


  void inithistos();
  void savehistos(char*);

  void Go(Long64_t numEvents);
  
 private: 
  char m_filename[256]; 
  int m_nev;
  TH2D *m_deltafrac;
  TH2D *m_deltalayerfrac;
  TH2D *m_deltalayerfraceng;
  TH2D *m_deltalayerfracengnorm;
  TH2D *m_hits;
  TH2D *m_layers;
  TH1D *m_delta;
  


  ClassDef(afterglow,1) 
    
};

#endif
