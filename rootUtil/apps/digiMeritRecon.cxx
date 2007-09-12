#include "rootUtil/CompositeEventList.h"
#include <TFile.h>
#include <TTree.h>
#include <TRandom.h>
#include "TSystem.h"


#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <vector>

#ifndef WIN32
#include <unistd.h>
#else
#include "apps/winutil/XGetopt.h"
#endif

void usage() {

  using std::cout;
  using std::cerr;
  using std::endl;

  std::string m_theApp("digiMeritRecon.exe");

  cout << endl
       << "digiMeritRecon.exe" << endl
       << "Makes a pointer skim for digi, merit and recon files" << endl
       << endl;
  
  cout << "Usage:" << endl
       << "\t" << m_theApp << " [options] [input] -o <output>" << endl 
       << endl
       << "\tINPUT" << endl
       << "\t   -r <reconFile>    : recon ROOT files" << endl
       << "\t   -d <digiFile>     : digi ROOT files" << endl
       << "\t   -S <svacFile>     : svac ROOT files" << endl
       << "\t   -m <meritFile>    : merit ROOT files" << endl 
       << endl
       << "\t   -o <output>       : output file" << endl
       << endl
       << "\tOPTIONS for all jobs" << endl
       << "\t   -h                : print this message" << endl
       << "\t   -n <nEvents>      : run over <nEvents>" << endl
       << "\t   -s <startEvent>   : start with event <startEvent>" << endl
       << "\t   -R <prescale>     : Apply random prescale" << endl    
       << endl;
}
  

int main(int argn, char** argc) {
#ifdef WIN32
    gSystem->Load("libTreePlayer.dll");
#endif


  using std::cout;
  using std::cerr;
  using std::endl;
  using std::string;

  std::vector<string> digiFiles; 
  std::vector<string> meritFiles;
  std::vector<string> reconFiles;
  std::vector<string> svacFiles; 

  string outFile;

  UInt_t optval_s(0);
  UInt_t optval_n(0);

  UInt_t prescale(0);

  // parse options
  char* endPtr;  
  int opt;
  while ( (opt = getopt(argn, argc, "ho:d:r:S:m:n:s:R:")) != EOF ) {
    switch (opt) {
    case 'h':   // help      
      usage();
      return 1;
    case 'o':   //  output
      outFile = string(optarg);
      break;
    case 'd':   // digi files
      digiFiles.push_back(string(optarg));
      break;
    case 'r':   // recon files
      reconFiles.push_back(string(optarg));
      break;
    case 'S':   // Svac files
      svacFiles.push_back(string(optarg));
      break;
    case 'm':   // Merit files
      meritFiles.push_back(string(optarg));
      break;
    case 'n':   // number of events
      optval_n = strtoul( optarg, &endPtr, 0 );
      break;
    case 's':   // start event
      optval_s = strtoul( optarg, &endPtr, 0 );
      break;
    case 'R':   // start event
      prescale = strtoul( optarg, &endPtr, 0 );
      break;
    case '?':
      usage();
      return 2;
    default:
      cerr << opt << " not parsable..." << endl;
      cerr << "Try " << "digiMeritRecon.exe" << " -h" << endl;
      return 2;      
    }
  }

  if ( argn - optind > 0 ) {
    cerr << "digiMeritRecon.exe" << " only takes options, not bare arguments" << endl
	 << "Try " << "digiMeritRecon.exe" << " -h" << endl;
    return 3;
  }

  // First file only
  std::vector<TTree*> trees;

  CompositeEventList p;

  Long64_t nEvt(0);
  if ( digiFiles.size() > 0 ) {
    TFile* f = TFile::Open(digiFiles[0].c_str(),"READ");
    TTree* t = static_cast<TTree*>(f->Get("Digi"));
    nEvt = t->GetEntries();
    t->SetBranchStatus("*",0);
    trees.push_back(t);
    p.addComponent("Digi_");
  }
  if ( reconFiles.size() > 0 ) {
    TFile* f = TFile::Open(reconFiles[0].c_str(),"READ");
    TTree* t = static_cast<TTree*>(f->Get("Recon"));
    nEvt = t->GetEntries();
    t->SetBranchStatus("*",0);
    trees.push_back(t);
    p.addComponent("Recon_");
  }
  if ( svacFiles.size() > 0 ) {
    TFile* f = TFile::Open(svacFiles[0].c_str(),"READ");
    TTree* t = static_cast<TTree*>(f->Get("Output"));
    nEvt = t->GetEntries();
    t->SetBranchStatus("*",0);
    trees.push_back(t);
    p.addComponent("Svac_");
  }
  if ( meritFiles.size() > 0 ) {
    TFile* f = TFile::Open(meritFiles[0].c_str(),"READ");
    TTree* t = static_cast<TTree*>(f->Get("MeritTuple"));
    nEvt = t->GetEntries();
    t->SetBranchStatus("*",0);
    trees.push_back(t);
    p.addComponent("Merit_");
  }
  
  if ( outFile.length() < 2 ) {
    usage();
    return 2;
  }
  TFile* fout = p.makeFile(outFile.c_str(),"RECREATE");
  TRandom r;

  for ( Long64_t iEvt(0); iEvt < nEvt; iEvt++ ) {
    Double_t accept = r.Rndm();
    if ( (accept * (Double_t)prescale) > 1.0 ) {
      continue;
    }
    for ( UInt_t iTree(0); iTree < trees.size(); iTree++ ) {
      TTree* t = trees[iTree];
      t->LoadTree(iEvt);      
    }
    p.fillEvent(trees);
  }
  p.fillMeta();
  fout->Write();

  return 0;
}
