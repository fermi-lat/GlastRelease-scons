/** 
 * @file ntupleWriterSvc_load.cxx
 *
 * @brief Implementation of <Package>_load routine.
 * This routine is needed for forcing the linker to load all the 
 * components of the library. 
 * 
 * $Header$
 */

#include "GaudiKernel/DeclareFactoryEntries.h"

DECLARE_FACTORY_ENTRIES(ntupleWriterSvc) {
    DECLARE_SERVICE(ntupleWriterSvc);
    DECLARE_ALGORITHM(WriteTupleAlg);
} 

