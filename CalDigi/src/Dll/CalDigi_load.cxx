/** 
* @file GlastDigi_load.cpp
* @brief This is needed for forcing the linker to load all components
* of the library.
*
*  $Header$
*/

#include "GaudiKernel/DeclareFactoryEntries.h"
#include "GaudiKernel/IToolFactory.h"

#define DLL_DECL_TOOL(x)       extern const IToolFactory& x##Factory; x##Factory.addRef();


DECLARE_FACTORY_ENTRIES(CalDigi) {
    DECLARE_ALGORITHM( CalDigiAlg );
    DLL_DECL_TOOL( LinearTaper   );
    DLL_DECL_TOOL( OnePlusExpTaper   );
    DLL_DECL_TOOL( CalDigiRandom   );
} 

