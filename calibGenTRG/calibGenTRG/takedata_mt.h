/** 
*@class takedata_mt
*
*@brief Root Tree Ananlysis takedata_mt as a root tree analysis
*
* This class is for Trigger Root Tree Ananlysis of the time structure of a run
*
* Version 0.1 24-August-2006 Martin Kocian Creation
*/

#ifndef TAKEDATA_MT_H
#define TAKEDATA_MT_H 1

#if !defined(__CINT__)
// Need these includes if we wish to compile this code
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TMath.h>
#include <TStyle.h>
#include <TH2.h>
#include <TProfile.h>
#include "calibGenTRG/RootTreeAnalysis.h"
#else  // for interactive use
#include <TH1D.h>
class RootTreeAnalysis;
#endif 

#include  <iostream> 
#include  <iomanip> 
#include  <fstream> 

class takedata_mt : public RootTreeAnalysis {

public :

  takedata_mt(); 

  virtual ~takedata_mt() {;}
  int readParameterFile(char *);
  void writeoutresults(char*);


  void inithistos();
  void savehistos(char*);

  void Go(Long64_t numEvents);
  
 private: 
  char m_filename[256]; 
  int m_gid;
  int m_nev;
  int m_tack_cal,m_tack_tkr,m_tack_acd;
  int m_evPerStep,m_stepNumber,m_temidlist[16];
  int m_totalcount[16],m_totalsquare[16];
  double m_totalgood[16];
  TH1D *m_tkrocc[16], *m_adchist[16], *m_acdadchist;
  
  string m_testnr; 
 


  ClassDef(takedata_mt,1) 
    
};

#endif
