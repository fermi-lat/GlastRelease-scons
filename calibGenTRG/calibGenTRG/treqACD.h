/** 
*@class treqACD
*
*@brief Root Tree Ananlysis treqACD as a root tree analysis
*
* This class is for Trigger Root Tree Ananlysis of the time structure of a run
*
* Version 0.1 24-August-2006 Martin Kocian Creation
*/

#ifndef TREQACD_H
#define TREQACD_H 1

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

class AcdId;
class TestReport;

class treqACD : public RootTreeAnalysis {

public :

  treqACD(); 

  virtual ~treqACD() {;}

  void inithistos();
  void setParameters(int tkrdelay,int acddelay);
  void writeoutresults(const char* reportname, const char* filename);
  void fit(TH1D* histo,double& mean, double &emean, double& sigma, double& esigma, int& nevents);
  void Go(Long64_t numEvents);
  void markError(char*, TestReport*);
  
 private: 
  char m_filename[256]; 
  int m_gid;
  double m_starttime;
  int m_nev;
  int m_status;
  int m_tkrdelay,m_acddelay;
  double m_mean,m_sigma;
  TH1D *m_allevents, *m_intersect, *m_singletile, *m_singletower, *m_singletowertile;
  TH1D *m_bytower[16];
  TH1D *m_bytile[108];
  
  string m_testnr; 
 


  ClassDef(treqACD,1) 
    
};

#endif
