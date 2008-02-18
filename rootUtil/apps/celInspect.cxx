
#include <rootUtil/CompositeEventList.h>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>

#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <vector>

#ifndef WIN32
#include <unistd.h>
#else
#include "apps/winutil/XGetopt.h"
#endif

using std::cout ;
using std::cerr ;
using std::endl ;
using std::string ;

void usage()
 {
  std::string m_theApp("celInspect.exe") ;

  std::cout << std::endl
       << m_theApp << std::endl
       << "Prints information about a Composite Event List" << std::endl
       << std::endl ;
  
  std::cout << "Usage:" << std::endl
       << "\t" << m_theApp << " [options] CelFile" << std::endl 
       << std::endl
       << std::endl
       << "\tOPTIONS for all jobs" << std::endl
       << "\t   -h                : print this message" << std::endl
       << std::endl ;
 }
  

int main(int argn, char** argc)
 {
# ifdef WIN32
    gSystem->Load("libTreePlayer.dll");
# endif

  std::string printFlags;
  UInt_t optval_s(0);
  UInt_t optval_n(10);

  // parse options
  char* endPtr;  
  int opt;
  while ( (opt = getopt(argn, argc, "hp:n:s:")) != EOF )
   {
    switch (opt)
     {
    case 'h':   // help      
      usage();
      return 1;
    case '?':
      usage();
      return 2;
    default:
      cerr << opt << " not parsable..." << endl;
      return 2;      
     }
   }

  if ( argn - optind == 0 )
   {
    cerr << "celInspect.exe" << " require an input file" << endl
         << "Try " << "celInspect.exe" << " -h" << endl ;
    return 3 ;
   }

  if ( argn - optind > 1 )
   {
    cerr << "celInspect.exe" << " accept only one input file" << endl
         << "Try " << "celInspect.exe" << " -h" << endl ;
    return 4 ;
   }

  const char * fileName = argc[optind] ;
  CompositeEventList cel(fileName) ;
  if ( cel.isOk()==kFALSE )
   {
    cerr << "Can't open file " << fileName << endl;
    return 5 ;
   }
  cel.printInfo() ;

  return 0 ;
 }
