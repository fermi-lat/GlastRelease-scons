// Tool to read a MOOT configuration and convert it to ROOT for later processing
//
//
//
//
//
// stl
#include <string>
#include <fstream>
#include <iostream>

#ifndef WIN32
#include <unistd.h>
#else
#include "facilities/XGetopt.h"
#endif
//

// local
#include "../src/fsw/fsw_cdb.h"

void usage() {
  std::cout << "printCDB.exe -i <CDM_NAME> -o <OUTFILE>" << std::endl
	    << std::endl
	    << "     CDM_NAME             : Name of CDM" << std::endl
	    << "     OUTFILE              : Name of output file" << std::endl
	    << std::endl;  
}

//
int main(int argn, char** argc) {
  
  // parse options
  // char* endPtr;  
  int opt;

  std::string outFileName;
  std::string cdm_name;

#ifdef WIN32
  while ( (opt = facilities::getopt(argn, argc, "hi:o:")) != EOF ) {
#else
  while ( (opt = getopt(argn, argc, "hi:o:")) != EOF ) {
#endif

    switch (opt) {
    case 'h':   // help      
      usage();
      return 1;
    case 'o':   //  output
      outFileName = std::string(optarg);
      break;
    case 'i':   //  cdm_name
      cdm_name = std::string(optarg);
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

  if ( outFileName.empty() ) {
    std::cerr << "No output file specified." << std::endl;
    usage();
    return 2;
  } else {
    std::cout << "Output file :  " << outFileName << std::endl;
  }
 
  if ( cdm_name.empty() ) {
    std::cerr << "No input CDM_NAME specified." << std::endl;
    usage();
    return 2;    
  } else {
    std::cout << "CDM_NAME    :  " <<  cdm_name << std::endl;
  }
  
  configData::fsw_cdb* cdb = configData::fsw_cdb::read_cdb(cdm_name.c_str(),0);
  if ( cdb == 0 ) {
    std::cerr << "Failed to read CDB " << cdm_name << std::endl;
    return 3;
  }
  
  if ( cdb->printToStream( std::cout ) != 0 ) {
    return 4;
  }

  if ( cdb->writeToXmlFile( outFileName.c_str() ) != 0 ) {
    std::cerr << "Failed to print CDB " << cdm_name 
	      << " to file " << outFileName << std::endl;
    return 4;
  }
  
  

  return 0;
}







