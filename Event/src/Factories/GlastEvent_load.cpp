
/*! \file GlastEvent_load.cpp
\brief based upon LHCbEvent_load.cpp by Markus Frank
*/

#define GLASTEVENT_FACTORIES_GLASTEVENT_LOAD_CPP  1

extern void   HitInstanciation();
extern void   IrfInstantiation();
extern void   MCInstantiation();
extern void   TopInstantiation();
extern void   RawInstantiation();

void GlastEventFactories_load()     {
    //HitInstanciation ();
    IrfInstantiation();
    MCInstantiation ();
    TopInstantiation ();
    RawInstantiation ();
}

extern "C" void GlastEventFactories_loadRef()  {
  GlastEventFactories_load();
}
