//====================================================================
// GlastEvent_load.cpp
//--------------------------------------------------------------------
//
//	Package   : GlastEvent
//
//	Author    : Markus Frank
//  History   :
// +---------+----------------------------------------------+---------
// |    Date |                 Comment                      | Who     
// +---------+----------------------------------------------+---------
// | 21/07/99| Initial version                              | MF
// +---------+----------------------------------------------+---------
//====================================================================
#define GLASTEVENT_FACTORIES_GLASTEVENT_LOAD_CPP  1

extern void   HitInstanciation();
extern void   MCInstanciation();
extern void   TopInstanciation();

void GlastEventFactories_load()     {
    //HitInstanciation ();
    MCInstanciation ();
    TopInstanciation ();
}

extern "C" void GlastEventFactories_loadRef()  {
  GlastEventFactories_load();
}
