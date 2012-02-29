
//
// stl
#include <string>
#include <iostream>
#ifndef WIN32
#include <unistd.h>
#else
#include "facilities/XGetopt.h"
using facilities::getopt;
#endif
//
// ROOT-io
#include "TFile.h"
#include "TTree.h"

// local
#include "configData/ConfigBase.h"
#include "configData/ConfigCompare.h"
#include "configData/ConfigReport.h"

void usage() {
  std::cout << "RootConfigComapre.exe [options] file1.root file2.root" << std::endl
	    << "\tcompares two configuration files" << std::endl
	    << "\t\tOptions:" << std::endl
	    << "\t-r\tPrint config report" << std::endl
	    << "\t-n\tNo comparison" << std::endl
	    << "\t-1\tOnly print one difference per register type" << std::endl
	    << "\t-f\tFull comparison, include threshold registers, implies -1" << std::endl
	    << "\t-o\tOutput file name stub" << std::endl;    
}

//
int main(int argn, char** argc) {
  
  // parse options
  //char* endPtr;  

  Bool_t onlyOne(kFALSE);
  Bool_t fullCompare(kFALSE);
  Bool_t report(kFALSE);
  Bool_t noDiff(kFALSE);
  std::string outputName;
  std::string outputPrefix;

  int opt;
  while ( (opt = getopt(argn, argc, "hf1rno:")) != EOF ) {
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
    case 'o':   //  output
      outputPrefix = std::string(optarg);
      break;
    default:
      std::cout << opt << " not parsable..." << std::endl;
      usage();
      return 2;  
    }
  }
  std::string fn1 = argc[optind];  
  outputName = outputPrefix+"report.root";

  TFile* f1 = TFile::Open(fn1.c_str());
  TTree* t1 = (TTree*)f1->Get("Config");

  ConfigBase c1(t1);
  c1.GetEntry(0); c1.LoadTree(0); 

  if ( ! noDiff ) {
    TList l;
    std::string fn2 = argc[optind+1];
    TFile* f2 = TFile::Open(fn2.c_str());
    TTree* t2 = (TTree*)f2->Get("Config");
    ConfigBase c2(t2);
    c2.GetEntry(0); c2.LoadTree(0); 
    ConfigCompare comp;
    comp.compare(c1, c2, l, fullCompare, onlyOne);
    TFile f(outputName.c_str(),"RECREATE");
    l.Write();
    f.Close();
  }

  if ( report ) {
    TList l;
    ConfigReport rep;
    rep.report(c1,l,"",std::cout);
    char* openOption = (! noDiff) ? "UPDATE" : "RECREATE";
    TFile f(outputName.c_str(),openOption);
    l.Write();
    f.Close();
  }

  return 0;
}







