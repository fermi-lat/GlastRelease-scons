
/*! \file IrfInstantiation.cpp
\brief based upon mcInstantiation.cpp by Markus Frank available within the LHCbEvent package

This file instanciates concretely the implementation of all of these
classes so that they may be included within factories in the DLL.

Original author: Sawyer Gillespie, hgillesp@u.washington.edu
*/

// Include files
#include "Gaudi/Kernel/ObjectFactory.h"
#include "GlastEvent/TopLevel/ObjectVector.h"
#include "GlastEvent/TopLevel/ObjectList.h"

#define _ImplementHitContainedFactories(x)  \
  _ImplementContainedObjectFactory( x )    \
  _ImplementDataObjectFactory( x##Vector ) \
  _ImplementDataObjectFactory( x##List )

#define DLL_DECL_CONTAINEDOBJECTFACTORY(x)  \
  DLL_DECL_OBJECTFACTORY( x ) \
  DLL_DECL_OBJECTFACTORY( x##Vector )\
  DLL_DECL_OBJECTFACTORY( x##List )

#include "GlastEvent/TopLevel/EventModel.h"
using namespace GlastEvent;

// ====================================================================
// Object factory implementation for objects of class IrfAcdHit 
// ====================================================================
#include "GlastEvent/Irf/IrfAcdHit.h"
_ImplementHitContainedFactories(IrfAcdHit)

// ====================================================================
// Object factory implementation for objects of class IrfCalHit
// ====================================================================
#include "GlastEvent/Irf/IrfCalHit.h"
_ImplementHitContainedFactories(IrfCalHit)


// ====================================================================
// Object factory implementation for objects of class IrfTkrHit
// ====================================================================
#include "GlastEvent/Irf/IrfTkrHit.h"
_ImplementHitContainedFactories(IrfTkrHit)

// ====================================================================
// Object factory implementation for objects of class IrfTkrLayer
// ====================================================================
#include "GlastEvent/Irf/IrfTkrLayer.h"
_ImplementHitContainedFactories(IrfTkrLayer)

void IrfInstantiation()  {
    DLL_DECL_CONTAINEDOBJECTFACTORY( IrfAcdHit );
    DLL_DECL_CONTAINEDOBJECTFACTORY( IrfCalHit );
    DLL_DECL_CONTAINEDOBJECTFACTORY( IrfTkrHit );
    DLL_DECL_CONTAINEDOBJECTFACTORY( IrfTkrLayer );
}
