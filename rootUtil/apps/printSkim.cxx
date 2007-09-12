#include "rootUtil/PointerSkim.h"
#include <TFile.h>
#include <TTree.h>
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

  std::string m_theApp("printSkim.exe");

  cout << endl
       << m_theApp << endl
       << "Prints information about a pointer skim file" << endl
       << endl;
  
  cout << "Usage:" << endl
       << "\t" << m_theApp << " [options] inputFiles" << endl 
       << endl
       << endl
       << "\tOPTIONS for all jobs" << endl
       << "\t   -h                : print this message" << endl
       << "\t   -p <flags>        : what to print"      << endl
       << "\t   -n <nEvents>      : run over <nEvents>" << endl
       << "\t   -s <startEvent>   : start with event <startEvent>" << endl
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


  std::string printFlags;
  UInt_t optval_s(0);
  UInt_t optval_n(10);

  // parse options
  char* endPtr;  
  int opt;
  while ( (opt = getopt(argn, argc, "hp:n:s:")) != EOF ) {
    switch (opt) {
    case 'h':   // help      
      usage();
      return 1;
    case 'p':   //  output
      printFlags = string(optarg);
      break;
    case 'n':   // number of events
      optval_n = strtoul( optarg, &endPtr, 0 );
      break;
    case 's':   // start event
      optval_s = strtoul( optarg, &endPtr, 0 );
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

  if ( argn - optind == 0 ) {
    cerr << "printSkim.exe" << " requires some input files" << endl
	 << "Try " << "printSkim.exe" << " -h" << endl;
    return 3;
  }

  PointerSkim p;

  for ( int idx = optind; idx < argn; idx++ ) {
    const char* fileName = argc[idx];
    cout << fileName << endl;
    TFile* f = p.openFile(fileName);
    if ( 0 == f ) {
      cerr << "Can't open file " << fileName << endl;
      return 4;
    }
    p.printout(printFlags.c_str(),optval_n,optval_s);
  }
  

  return 0;
}
