// 
//  Original author: Sawyer Gillespie
//                   hgillesp@u.washington.edu
//

// *********
//  This file instanciates concretely the implementation of all of these
//  classes so that they may be included within factories in the DLL.
//

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
// Object factory implementation for objects of class MCACDHit 
// ====================================================================
#include "GlastEvent/MonteCarlo/MCACDHit.h"
_ImplementHitContainedFactories(MCACDHit)

// ====================================================================
// Object factory implementation for objects of class MCCalorimeterHit
// ====================================================================
#include "GlastEvent/MonteCarlo/MCCalorimeterHit.h"
_ImplementHitContainedFactories(MCCalorimeterHit)


void MCInstanciation()  {
    DLL_DECL_CONTAINEDOBJECTFACTORY( MCACDHit );
    DLL_DECL_CONTAINEDOBJECTFACTORY( MCCalorimeterHit );
}
