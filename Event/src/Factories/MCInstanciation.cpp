
/*! \file MCInstanciation.cpp
\brief based upon mcInstanciation.cpp by Markus Frank available within the LHCbEvent package

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
// Object factory implementation for objects of class McParticle 
// ====================================================================
#include "GlastEvent/MonteCarlo/McParticle.h"
_ImplementHitContainedFactories(McParticle)

// ====================================================================
// Object factory implementation for objects of class McPositionHit
// ====================================================================
#include "GlastEvent/MonteCarlo/McPositionHit.h"
_ImplementHitContainedFactories(McPositionHit)


// ====================================================================
// Object factory implementation for objects of class McVertex
// ====================================================================
#include "GlastEvent/MonteCarlo/McVertex.h"
_ImplementHitContainedFactories(McVertex)

// ====================================================================
// Object factory implementation for objects of class McIntegratingHit
// ====================================================================
#include "GlastEvent/MonteCarlo/McIntegratingHit.h"
_ImplementHitContainedFactories(McIntegratingHit)

// ====================================================================
// Object factory implementation for objects of class CsIData 
// ====================================================================
void MCInstantiation()  {
    DLL_DECL_CONTAINEDOBJECTFACTORY( McParticle );
    DLL_DECL_CONTAINEDOBJECTFACTORY( McPositionHit );
    DLL_DECL_CONTAINEDOBJECTFACTORY( McVertex );
    DLL_DECL_CONTAINEDOBJECTFACTORY( McIntegratingHit );
}
