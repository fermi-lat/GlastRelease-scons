/** 

*  $Header$
*/

#include "GaudiKernel/DeclareFactoryEntries.h"
#include "GaudiKernel/IToolFactory.h"

#define DLL_DECL_TOOL(x) extern const IToolFactory& x##Factory; x##Factory.addRef();


DECLARE_FACTORY_ENTRIES(HepRepCorba) {
  DLL_DECL_TOOL(RegisterCorba);
}

#include "GaudiKernel/LoadFactoryEntries.h"

LOAD_FACTORY_ENTRIES(HepRepCorba)
