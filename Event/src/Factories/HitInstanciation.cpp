// $Id$
// 
//  Original author: Sawyer Gillespie
//                   hgillesp@u.washington.edu
//

// *********
//  This file instanciates concretely the Gsim implementation of all of these
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
// Object factory implementation for objects of class GsimACDTileHits
// ====================================================================
#include "GlastEvent/Hits/ACDhit.h"
_ImplementHitContainedFactories(ACDhit)

#include "GlastEvent/Hits/MCCalorimeterHit.h"
_ImplementHitContainedFactories(MCCalorimeterHit)

// ====================================================================
// Object factory implementation for objects of class GsimCalorimeterHits
// ====================================================================
//#include "GlastEvent/GsimImp/GsimCalorimeterHits.h"
//_ImplementDataObjectFactory(GsimCalorimeterHits)

// ====================================================================
// Object factory implementation for objects of class GsimCalorimeterLogHits
// ====================================================================
//#include "GlastEvent/GsimImp/GsimCalorimeterLogHits.h"
//_ImplementGsimImpContainedFactories(GsimCalorimeterLogHits)

// ====================================================================
// Object factory implementation for objects of class GsimGlastHits
// ====================================================================
//#include "GlastEvent/GsimImp/GsimGlastHits.h"
//_ImplementDataObjectFactory(GsimGlastHits)

// ====================================================================
// Object factory implementation for objects of class GsimSiLayerHits
// ====================================================================
//#include "GlastEvent/GsimImp/GsimSiLayerHits.h"
//_ImplementDataObjectFactory(GsimSiLayerHits)

// ====================================================================
// Object factory implementation for objects of class GsimSiStripHits
// ====================================================================
//#include "GlastEvent/GsimImp/GsimSiStripHits.h"
//_ImplementGsimImpContainedFactories(GsimSiStripHits)

// ====================================================================
// Object factory implementation for objects of class GsimTowerHits
// ====================================================================
//#include "GlastEvent/GsimImp/GsimTowerHits.h"
//_ImplementDataObjectFactory(GsimTowerHits)

// ====================================================================
// Object factory implementation for objects of class GsimTrackerHits
// ====================================================================
//#include "GlastEvent/GsimImp/GsimTrackerHits.h"
//_ImplementDataObjectFactory(GsimTrackerHits)

// ====================================================================
// Object factory implementation for objects of class MCHitBase
// Objects cannot be instantiated!
// ===================================================================
//_ImplementDataObjectFactory( GsimGlastHits )
//_ImplementDataObjectFactory( GsimTowerHits )
//_ImplementDataObjectFactory( GsimTrackerHits )
//_ImplementDataObjectFactory( GsimCalorimeterHits )
//_ImplementDataObjectFactory( GsimSiLayerHits )

void HitInstanciation()  {
  DLL_DECL_CONTAINEDOBJECTFACTORY( ACDhit );
  DLL_DECL_CONTAINEDOBJECTFACTORY( MCCalorimeterHit );
  //DLL_DECL_CONTAINEDOBJECTFACTORY( GsimCalorimeterLogHits );
  //DLL_DECL_CONTAINEDOBJECTFACTORY( GsimSiStripHits );
}
