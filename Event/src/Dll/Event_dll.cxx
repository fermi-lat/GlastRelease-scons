
#include "GaudiKernel/LoadFactoryEntries.h"

LOAD_FACTORY_ENTRIES(Event)


#define define_dll_export 

#include "Event/TopLevel/EventModel.h"

//#include "GaudiKernel/ClassID.h"

//#define dllexport __declspec( dllexport )

//dllexport const CLID& CLID_Event = 110;