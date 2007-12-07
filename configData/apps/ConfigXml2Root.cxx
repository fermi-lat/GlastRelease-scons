
//
// stl
#include <string>
#include <fstream>

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
#include "configData/LatConfig.h"
#include "configData/ConfigReader.h"

void usage() {
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

#ifdef WIN32
  while ( (opt = facilities::getopt(argn, argc, "ho:p:a:b")) != EOF ) {
#else
  while ( (opt = getopt(argn, argc, "ho:p:a:b")) != EOF ) {
#endif
    switch (opt) {
    case 'h':   // help      
      usage();
      return 1;
    case 'o':   //  output
      outputPrefix = std::string(optarg);
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
    }else{
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
    reader.readTopLvl(itrIO->first.c_str(),paramDir,bcast);
    config.latch();
    
    t->Write();
    fOut.Close();
  }

  return 0;
}







