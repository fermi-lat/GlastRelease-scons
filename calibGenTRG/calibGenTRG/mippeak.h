/** 
*@class mippeak
*
*@brief Root Tree Ananlysis mippeak as a root tree analysis
*
* This class is for Trigger Root Tree Ananlysis of the time structure of a run
*
* Version 0.1 24-August-2006 Martin Kocian Creation
*/

#ifndef MIPPEAK_H
#define MIPPEAK_H 1

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

class mippeak : public RootTreeAnalysis {

public :

  mippeak(); 

  virtual ~mippeak() {;}
  void writeoutresults(char*);


  void inithistos();
  void savehistos(char*);

  void Go(Long64_t numEvents);
  
 private: 
  char m_filename[256]; 
  int m_gid;
  double m_starttime;
  int m_nev;
  int m_evPerStep,m_stepNumber,m_temidlist[16];
  double m_mp,m_mp_err;
  double m_lanwid, m_lanwid_err;
  double m_gaussig, m_gaussig_err;
  double m_chi2;
  TH1D *m_acdadchist;
  
  string m_testnr; 
 


  ClassDef(mippeak,1) 
    
};

#endif
