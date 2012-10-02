// $Header$

// Include files
#include "GaudiKernel/SmartIF.h"
#include "GaudiKernel/Bootstrap.h"
#include "GaudiKernel/IAppMgrUI.h"
#include "GaudiKernel/IProperty.h"
#include "GaudiKernel/IJobOptionsSvc.h"
#include "GaudiKernel/System.h"
// system stuff
#ifdef WIN32
# include <direct.h>
#else
# include <unistd.h>
#endif
#include <time.h>
#include <stdio.h>

#include "facilities/commonUtilities.h"
#include "facilities/Util.h"

//------------------------------------------------------------------------------
//
//  Package    : GlastPolicy
//
//  Description: generic Gaudi Main Program
//
//------------------------------------------------------------------------------
/**
Standard Glast-Gaudi main.
Specify the job options file location by:
1) specification, local or global, on command line.
2) content of env var JOBOPTIONS
3) src/jobOptions

*/
// declare function that may reduce the execution priority
void setPriority();
#ifdef WIN32
#define chdir _chdir
#endif

void current_time(std::ostream& out=std::cout)
{   
        static bool first=true;
        static time_t start;
        if(first){ first=false; ::time(&start);}
        time_t aclock;
        ::time( &aclock );   
        char tbuf[25]; ::strncpy(tbuf, asctime( localtime( &aclock ) ),24);
        tbuf[24]=0;
        out<<  "Current time: " << tbuf
            << " ( "<< ::difftime( aclock, start) <<" s elapsed)" << std::endl;
}


int main( int argn, char** argc) {

  facilities::commonUtilities::setupEnvironment();
    
    //    std::string joboptions_file="src/jobOptions.txt"; // default
  std::string joboptions_file;
#ifndef SCons    
  joboptions_file="src/jobOptions.txt"; // default
#else
#ifdef PACKAGE_NAME
  std::string pkgName(PACKAGE_NAME); 
  joboptions_file = facilities::commonUtilities::getJobOptionsPath(pkgName) +
    "/jobOptions.txt";
#else
  std::cerr << "SCons compile did not specify package name! ";
  return 1;
#endif
#endif
    
  // check for env var
  std::string jobstring = facilities::commonUtilities::getEnvironment("JOBOPTIONS");

  // Priority goes to command arg, but it must be non-empty
  // (If called from wrapper script may see empty first argument)
  if( argn>1 ) {
    std::string arg1(argc[1]);
    if (arg1.size() > 0) {
      joboptions_file = argc[1]; // priority to command arg.
    }
  }
  else if(jobstring.size() > 0 ) { joboptions_file = jobstring; }
 

  // translate env variables if any
  try {
    facilities::Util::expandEnvVar(&joboptions_file);
  }
  catch (facilities::Untranslatable ex) {
    std::cerr << "Unable to expand job options name " << joboptions_file 
              << std::endl;
    return 1;
  }

  std::cerr << "Starting Glast-Gaudi job with job options file " 
            << joboptions_file << std::endl;
  current_time(std::cerr);
  const char* newdir = ::getenv("GLEAM_CHDIR");
  if( newdir !=0) {
    std::cerr << "Changing default directory to: " << newdir << std::endl;
    ::chdir(newdir);
  }
  
  
  // Create an instance of an application manager
  IInterface* iface = Gaudi::createApplicationMgr();
  
  SmartIF<IProperty>     propMgr ( iface );
  SmartIF<IAppMgrUI>     appMgr  ( iface );
  
  if( !appMgr.isValid() || !propMgr.isValid() ) {
    std::cout << "Fatal error while creating the ApplicationMgr " << std::endl;
    return 1;
  }
  
  
  // Set properties of algorithms and services
  propMgr->setProperty( "JobOptionsPath", joboptions_file );
  
  // Run the application manager and process events
#ifdef WIN32
  setPriority();
#endif
  StatusCode status = appMgr->run();

    // All done - exit with 0 if success.
    iface->release();
    current_time(std::cerr);
    return (status.isFailure()? 1 : 0);
    
}
