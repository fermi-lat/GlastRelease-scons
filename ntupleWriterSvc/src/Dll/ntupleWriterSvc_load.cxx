/** 
 * @file ntupleWriterSvc_load.cxx
 *
 * @brief Implementation of <Package>_load routine.
 * This routine is needed for forcing the linker to load all the 
 * components of the library. 
 * 
 * @author H. Gillespie
 * 
 * $Header$
 */

#include "GaudiKernel/ICnvFactory.h"
#include "GaudiKernel/ISvcFactory.h"
#include "GaudiKernel/IAlgFactory.h"


#define DLL_DECL_SERVICE(x)    extern const ISvcFactory& x##Factory; x##Factory.addRef();
#define DLL_DECL_CONVERTER(x)  extern const ICnvFactory& x##Factory; x##Factory.addRef();
#define DLL_DECL_ALGORITHM(x)  extern const IAlgFactory& x##Factory; x##Factory.addRef();
#define DLL_DECL_OBJECT(x)     extern const IFactory& x##Factory; x##Factory.addRef();

//! Load all  services: 
void ntupleWriterSvc_load() {
    DLL_DECL_SERVICE( ntupleWriterSvc );
    DLL_DECL_ALGORITHM( WriteTupleAlg );
} 

extern "C" void ntupleWriterSvc_loadRef()    {
  ntupleWriterSvc_load();
}

