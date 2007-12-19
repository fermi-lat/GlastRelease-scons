
//
// stl
#include <string>
#include <iostream>
#ifndef WIN32
#include <unistd.h>
#else
#include "facilities/XGetopt.h"
#endif
//
// ROOT-io
#include "TFile.h"
#include "TTree.h"

// local
#include "configData/ConfigCompare.h"

void usage() {
  std::cout << "RootConfigComapre.exe [options] file1.root file2.root" << std::endl
	    << "\tcompares two configuration files" << std::endl
	    << "\t\tOptions:" << std::endl
	    << "\t-r\tPrint config report" << std::endl
	    << "\t-1\tOnly print one difference per register type" << std::endl
	    << "\t-f\tFull comparison, include threshold registers, implies -1" << std::endl;    
}

//
int main(int argn, char** argc) {
  
  // parse options
  //char* endPtr;  

  Bool_t onlyOne(kFALSE);
  Bool_t fullCompare(kFALSE);
  Bool_t report(kFALSE);
  Bool_t noDiff(kFALSE);

  int opt;
#ifdef WIN32
  while ( (opt = facilities::getopt(argn, argc, "hf1rn")) != EOF ) {
#else
  while ( (opt = getopt(argn, argc, "hf1rn")) != EOF ) {
#endif
    switch (opt) {
    case 'h':   // help      
      usage();
      return 1;
    case '1':
      onlyOne = kTRUE;
      break;
    case 'r':
      report = kTRUE;
      break;
    case 'f':
      fullCompare = kTRUE;
      onlyOne = kTRUE;
      break;
   case 'n':
      noDiff = kTRUE;
      break;
    default:
      std::cout << opt << " not parsable..." << std::endl;
      usage();
      return 2;  
    }
  }
  std::string fn1 = argc[optind];  

  TFile* f1 = TFile::Open(fn1.c_str());
  TTree* t1 = (TTree*)f1->Get("Config");

  ConfigCompare c1(t1);
  
  if ( ! noDiff ) {
    std::string fn2 = argc[optind+1];
    TFile* f2 = TFile::Open(fn2.c_str());
    TTree* t2 = (TTree*)f2->Get("Config");
    ConfigCompare c2(t2);
    c1.compare(c2,0,fullCompare,onlyOne);
  }

  if ( report ) {
    TList l;
    c1.report(l);
    TFile f("report.root","RECREATE");
    l.Write();
    f.Close();
  }

  return 0;
}







