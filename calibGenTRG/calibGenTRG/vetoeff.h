/** 
*@class vetoeff
*
*@brief Root Tree Analysis vetoeff as a root tree analysis
*
* This class is for analyzing veto efficiency runs.
*
* Version 0.1 17-October-2006 Martin Kocian Creation
*/

#ifndef VETOEFF_MT_H
#define VETOEFF_MT_H 1

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
#include <TProfile.h>
class RootTreeAnalysis;
#endif 

#include  <iostream> 
#include  <iomanip> 
#include  <fstream> 
#include <vector>

class vetoeff : public RootTreeAnalysis {

public :

  vetoeff(); 

  virtual ~vetoeff() {;}
  void readParameterFile(char *);
  void writeoutresults(char*);
  void setThreshold(float t){m_thresh=t;}
 

  void inithistos();
  void savehistos(char*);
  void Go(Long64_t numEvents);

  float efficiency(int,int);
  float efferror(int,int); 
  
 private: 
  char m_filename[256]; 
  int m_gid;
  int m_nev;
  int m_roicount,m_roicount1,m_roicount3;
  int m_roicount32,m_roicount35,m_roicountveto;
  int m_status;
  int m_32shadow,m_32nofilter,m_32tkrintersec,m_32noise;
  int m_baddigievents;
  int m_diginogem,m_gemdigi,m_gemnodigi;
  float m_thresh;
  std::vector<int> m_roi[16];
  int m_goodcount[16],m_badcount[16],m_goodcountev[16];
  int m_baddigicount[16];
  TH1D *m_goodtiles[16], *m_badtiles[16], *roicond;
  TH1D *m_baddigis[16];
  TH1D *digicond, *digiglt,*glt3;
  TH1D *m_tilet, *m_tiled;
  
 


  ClassDef(vetoeff,1) 
    
};

#endif
