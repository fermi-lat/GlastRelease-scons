#include "rootUtil/CelUtil.h"

#include <TFile.h>
#include <TTree.h>
#include <TList.h>
#include <TROOT.h>
#include "TSystem.h"


#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <assert.h>

#ifndef WIN32
#include <unistd.h>
#else
#include "facilities/XGetopt.h"
using facilities::getopt;
#endif

void usage() {


  std::string theApp("mergeSkim.exe");

  std::cout << std::endl
       << theApp << std::endl
       << "Merges composite event lists" << std::endl
       << std::endl;
  
  std::cout << "Usage:" << std::endl
	    << "\t" << theApp << "-o <output> [options] <inputFiles>" << std::endl 
	    << "\tOptions" << std::endl
	    << "\t\t-f\tForce output (overwrites file)"
	    << std::endl
	    << std::endl;
}
  

int main(int argn, char** argc) {
#ifdef WIN32
    gSystem->Load("libTreePlayer.dll");
#endif


  std::vector<std::string> inputFiles;
  std::string outFile;

  std::string theApp("mergeSkim.exe");
  
  bool optval_f(false);

  // parse options
  //char* endPtr;  
  int opt;
  while ( (opt = getopt(argn, argc, "ho:f")) != EOF ) {
    switch (opt) {
    case 'h':   // help      
      usage();
      return 1;
    case 'o':   // output
      outFile = std::string(optarg);
      break;
    case 'f':   // force output
      optval_f = true;
      break;
    case '?':
      usage();
      return 2;
    default:
      std::cerr << opt << " not parsable..." << std::endl;
      std::cerr << "Try " << theApp << " -h" << std::endl;
      return 2;      
    }
  }

  if ( outFile.length() == 0 ) {
    std::cerr << "No output file given" << std::endl
	      << "Try " << theApp << " -h" << std::endl;
    return 3;
  }
  if ( argn - optind == 0 ) {
    std::cerr << "No input files given" << std::endl
	      << "Try " << theApp << " -h" << std::endl;
    return 4;
  }

  // Init ROOT
  
  if ( ! gROOT->Initialized() ) {
    assert(0);
  }
  TList inputFileList;
  for ( int idx = optind; idx < argn; idx++ ) {
    const char* fileName = argc[idx];
    std::cout << fileName << std::endl;
    TFile* f = TFile::Open(fileName);
    if ( 0 == f ) {
      std::cerr << "Can't open file " << fileName << std::endl;
      return 5;
    }
    inputFileList.AddLast(f);
  }
  
  CompositeEventList* s = optval_f ? 
    CelUtil::mergeCelFiles(inputFileList,outFile.c_str(),"RECREATE") :
    CelUtil::mergeCelFiles(inputFileList,outFile.c_str());

  if ( 0 == s ) return 6;

  return 0;
}
