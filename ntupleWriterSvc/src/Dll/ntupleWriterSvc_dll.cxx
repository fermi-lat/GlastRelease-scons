/** 
 * @file ntupleWriterSvc_dll.cxx
 *
 * @brief Implementation of the DllMain routine.  
 * The DLL initialization must be done seperately for each DLL. 
 * 
 * $Header$
 */

#include "GaudiKernel/LoadFactoryEntries.h"
LOAD_FACTORY_ENTRIES(ntupleWriterSvc)


//void FATAL(const char * msg) {
//    std::cerr << "FATAL error from ntupleWriterSvc DLL: " << msg << std::endl;/
//}

//void WARNING(const char * t){ std::cerr << "WARNING from ntupleWriterSvc DLL: " << t << std::endl; }
