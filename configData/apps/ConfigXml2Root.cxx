// Tool to read a MOOT configuration and convert it to ROOT for later processing
//
//
//
//
//
// stl
#include <string>
#include <fstream>
#include <cstdlib>

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
#include "configData/LatConfig.h"
#include "configData/ConfigReader.h"

void usage() {
  std::cout << "ConfigXml2Root.exe [options] key1 ..." << std::endl
            << "\tgenerated a ROOT tree containing MOOT configuration information" << std::endl
            << "\t\tOptions"               << std::endl 
            << "\t-h\thelp (this message)" << std::endl
            << "\t-?\thelp (this message)" << std::endl
            << "\t-o\toutput file prefix"  << std::endl
            << "\t-M\tKey is for MOOT"     << std::endl
            << "\t-p\tparameter directory" << std::endl
            << "\t-a\tapp directory"       << std::endl
            << "\t-b\tbroadcast"           << std::endl;
}

//
int main(int argn, char** argc) {
  
  // parse options
  // char* endPtr;  
  int opt;

  std::string outputPrefix;
  std::string inputFile;
  std::string paramDir = "params";
  std::string appDir = "params";
  std::list<std::pair<std::string,std::string> > inAndOut;
  bool bcast(false);
  bool useMootKey=false;

  while ( (opt = getopt(argn, argc, "ho:Mp:a:b")) != EOF ) {
    switch (opt) {
    case 'h':   // help      
      usage();
      return 1;
    case 'o':   //  output
      outputPrefix = std::string(optarg);
      break;
    case 'M':   // mootKey
      useMootKey=true;
      break;
    case 'p':   // 
      paramDir = std::string(optarg);
      break;
    case 'a':   // 
      appDir = std::string(optarg);
      break;
    case 'b':
      bcast = true;
      break;
    case '?':
      usage();
      return 2;
    default:
      std::cerr << opt << " not parsable..." << std::endl;
      usage();
      return 2;  
    }
  }
  for ( int idx = optind; idx != argn; idx++ ) {
    char *endptr;
    strtoul(argc[idx],&endptr,0);
    std::string outName;
    std::string inName;
    if (endptr[0]=='\0'){ // user typed in a latc key
      inName=argc[idx];
      outName=outputPrefix+"MootConfig_"+inName+".root";
    }
    else {
      std::string baseName = argc[idx];
      if ( baseName.find(appDir) == 0 ) {
	baseName.erase(0,appDir.size()+1);
      }
      std::string outBase = baseName; outBase.replace(outBase.find(".xml"),5,std::string(".root"));
      inName = appDir; inName += "/"; inName += baseName;
      outName = outputPrefix; outName += outBase;
    }
    inAndOut.push_back(std::pair<std::string,std::string>(inName,outName));
  }
  

  const std::string nullString;
  LatConfig config("Config");  
  ConfigReader reader(config.tuple(),config.gem());
  
  for ( std::list<std::pair<std::string,std::string> >::const_iterator itrIO = inAndOut.begin();
	itrIO != inAndOut.end(); itrIO++ ) {

    TFile fOut(itrIO->second.c_str(),"RECREATE");
    TTree* t = config.makeTree(nullString);  
    reader.readTopLvl(itrIO->first.c_str(),paramDir,bcast,useMootKey);
    config.latch();
    
    t->Write();
    fOut.Close();
  }

  return 0;
}







