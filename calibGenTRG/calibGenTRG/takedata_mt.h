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
  int setParameters(int,int,int,int);
  void writeoutresults(char*);
  void initCalTuple(const char* filename);
  void initMeritTuple(const char* filename);
  void useKalmanCut(bool sw){m_useKalman=sw;}
  void useToTCut(bool sw){m_useToT=sw;}
  void useMipCut(bool sw){m_useMip=sw;}
  void useOneTrack(bool sw){m_useOneTrack=sw;}
  void setKalmanCut(float lower){m_kalmanLower=lower;}
  void setTotCut(float lower, float upper){m_ToTLower=lower;m_ToTUpper=upper;}
  void setMipLower(int n){m_MipLower=n;}
  void setPtMagLatCut(float lower, float upper){m_ptmaglatlower=lower;m_ptmaglatupper=upper;}
  void usePtMagLatCut(bool sw){m_useptmaglat=sw;}

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
  int m_MipLower;
  double m_totalgood[16];
  float m_calrange[16][8][12][2][4];
  float m_ptmaglat;
  TFile *m_caltuplefile;
  TTree* m_caltuple;
  TFile *m_merittuplefile;
  TTree* m_merittuple;
  TH1D *m_tkrocc[16], *m_adchist[16], *m_calratio[16], *m_acdadchist;
  bool m_useKalman;
  bool m_useToT;
  bool m_useMip;
  bool m_useOneTrack;
  float m_kalmanLower;
  float m_ToTLower;
  float m_ToTUpper;
  bool m_useptmaglat;
  float m_ptmaglatlower;
  float m_ptmaglatupper;
  string m_testnr; 
 


  ClassDef(takedata_mt,1) 
    
};

#endif
