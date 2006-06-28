/** 
* @file AncillaryDataRecon_load.cxx
* @brief This is needed for forcing the linker to load all components
* of the library.
*
*/

#include "GaudiKernel/DeclareFactoryEntries.h"

DECLARE_FACTORY_ENTRIES(AncillaryDataRecon) 
{
    DECLARE_ALGORITHM( AncillaryDataReconAlg );
}
