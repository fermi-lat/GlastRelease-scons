//====================================================================
//  FluxSvc_load.cpp
//--------------------------------------------------------------------
//
//  Package    : FluxSvc
//
//  Description: Implementation of <Package>_load routine. This routine 
//               is needed for forcing the linker to load all the components
//               of the library.. 
//
//====================================================================
#include "GaudiKernel/ICnvFactory.h"
#include "GaudiKernel/ISvcFactory.h"
#include "GaudiKernel/IAlgFactory.h"

#include "flux/ISpectrumFactory.h"

#define DLL_DECL_SERVICE(x)    extern const ISvcFactory& x##Factory; x##Factory.addRef();
#define DLL_DECL_CONVERTER(x)  extern const ICnvFactory& x##Factory; x##Factory.addRef();
#define DLL_DECL_ALGORITHM(x)  extern const IAlgFactory& x##Factory; x##Factory.addRef();

#define DLL_DECL_SPECTRUM(x)   extern const ISpectrumFactory& x##Factory; x##Factory.addRef();

void FluxSvc_load() {
  DLL_DECL_SERVICE( FluxSvc );

  // these are the spectra that we want to make available
  DLL_DECL_SPECTRUM( CHIMESpectrum);
  DLL_DECL_SPECTRUM( AlbedoPSpectrum);
  DLL_DECL_SPECTRUM( HeSpectrum);
  DLL_DECL_SPECTRUM( ProtonSpectrum);
  DLL_DECL_SPECTRUM( TrappedProtonSpectrum);
  DLL_DECL_SPECTRUM( GalElSpectrum);

  DLL_DECL_SPECTRUM( CrElectron);
  DLL_DECL_SPECTRUM( CrProton);
}

extern "C" void FluxSvc_loadRef() {
  FluxSvc_load();
}
