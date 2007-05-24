
//
// stl
#include <string>
#include <iostream>
#include <unistd.h>

//
// ROOT-io
#include "TFile.h"
#include "TTree.h"

// local
#include "configData/ConfigCompare.h"

void usage() {
}

//
int main(int argn, char** argc) {
  
  // parse options
  //char* endPtr;  
  int opt;

  while ( (opt = getopt(argn, argc, "h")) != EOF ) {
    switch (opt) {
    case 'h':   // help      
      usage();
      return 1;
    default:
      std::cout << opt << " not parsable..." << std::endl;
      usage();
      return 2;  
    }
  }
  std::string fn1 = argc[optind];  
  std::string fn2 = argc[optind+1];

  TFile* f1 = TFile::Open(fn1.c_str());
  TFile* f2 = TFile::Open(fn2.c_str());
  TTree* t1 = (TTree*)f1->Get("Config");
  TTree* t2 = (TTree*)f2->Get("Config");

  ConfigCompare c1(t1);
  ConfigCompare c2(t2);
  
  c1.compare(c2);

  return 0;
}







