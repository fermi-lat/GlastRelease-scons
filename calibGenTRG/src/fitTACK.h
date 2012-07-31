#ifndef fitTACK_h
#define fitTACK_h

#include <TH1D.h>
#include <TF1.h>
#include <fstream>
#include <iostream>
#include <TStyle.h>
#include <TCanvas.h>
#include <TFile.h>
#include <time.h>
#include "LangauFun.h"
#include "egfun.h"

#define MAXCAL 5
#define MAXTKR 0
#define MAXACD 5
#define STDACD 24
#define STDCAL 45
#define STDTKR 0
using namespace std;

class fittack{
public:
  fittack();
  virtual ~fittack(){}
  void settextinput(char*);
  void setrootoutput(char*);
  void setratiofilename(const char*);
  void setfirstlastpoint(int,int);
  void readtextfiles();
  void readrootfiles();
  TH1D* calhist(int i, int j){return adchist[i][j];}
  TH1D* acdhist(int i){return acdadchist[i];}
  TH1D* getAcdWaveform(){return acdwaveform;}
  void fit();
  int validate();
  void writereport(int status);
private:
  string textname;
  string rootoutput;
  int first,last,npoints;
  int evperstep[30],nev[30],tack_cal[30],tack_tkr[30],tack_acd[30];
  char gid[30][128];
  double stepacd,stepcal,steptkr;
  int totalcount[30][16],totalsquare[30][16];
  double totalgood[30][16];
  int totalcountsum[30],totalsquaresum[30];
  double totalgoodsum[30];
  string filename[30];
  string ratiofilename[30];
  bool useratiofile;
  double calfitmeansum,ecalfitmeansum;
  double calhifitmeansum,ecalhifitmeansum;
  double tkrfitmeansum,etkrfitmeansum;
  double calfitmean[16],ecalfitmean[16];
  double calhifitmean[16],ecalhifitmean[16];
  double tkrfitmean[16],etkrfitmean[16];
  double acdfitmean,eacdfitmean;
  int acdstatus,calstatus[16],tkrstatus[16];
  int calstatussum,tkrstatussum;
  TH1D *tkrocc[30][16], *adchist[30][16], *acdadchist[30], *calratio[30][16];
  TH1D  *adchistsum[30], *calratiosum[30];
  TH1D *acdwaveform, *calwaveform[16], *calhiwaveform[16], *tkrwaveform[30],*calratiowaveform[30];
  TH1D *calsumwaveform, *calhisumwaveform, *tkrsumwaveform, *calratiosumwaveform;
  ClassDef(fittack,1) 

};
#endif
