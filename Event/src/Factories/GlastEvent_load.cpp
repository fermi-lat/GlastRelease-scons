
/*! \file GlastEvent_load.cpp
\brief based upon LHCbEvent_load.cpp by Markus Frank
*/

#define GLASTEVENT_FACTORIES_GLASTEVENT_LOAD_CPP  1

extern void   HitInstanciation();
extern void   MCInstanciation();
extern void   TopInstantiation();
//extern void   RawInstanciation();

void GlastEventFactories_load()     {
    //HitInstanciation ();
    MCInstanciation ();
    TopInstantiation ();
}

extern "C" void GlastEventFactories_loadRef()  {
  GlastEventFactories_load();
}
