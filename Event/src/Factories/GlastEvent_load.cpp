
/*! \file GlastEvent_load.cpp
\brief based upon LHCbEvent_load.cpp by Markus Frank
*/

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
