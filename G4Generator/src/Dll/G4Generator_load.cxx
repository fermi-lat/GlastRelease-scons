//====================================================================
//====================================================================

#include "GaudiKernel/ICnvFactory.h"
#include "GaudiKernel/ISvcFactory.h"
#include "GaudiKernel/IAlgFactory.h"


#define DLL_DECL_SERVICE(x)    extern const ISvcFactory& x##Factory; x##Factory.addRef();
#define DLL_DECL_CONVERTER(x)  extern const ICnvFactory& x##Factory; x##Factory.addRef();
#define DLL_DECL_ALGORITHM(x)  extern const IAlgFactory& x##Factory; x##Factory.addRef();
#define DLL_DECL_OBJECT(x)     extern const IFactory& x##Factory; x##Factory.addRef();

//! Load all  services: 
void G4Generator_load() {
    DLL_DECL_ALGORITHM( G4Generator );

} 

extern "C" void G4Generator_loadRef()    {
  G4Generator_load();
}
// these are needed by old geometry routines: remove when fixed.
#include <iostream>
void WARNING(const char * msg){ std::cout  << msg << std::endl;}
void FATAL (const char* msg) {std::cerr << msg << std::endl; exit(-1);}