/*! \file DigiInstantiation.cpp
\brief based upon mcInstantiation.cpp by Markus Frank available within the LHCbEvent package

This file instanciates concretely the implementation of all of these
classes so that they may be included within factories in the DLL.

Original author: Heather Kelly
*/

// Include files
#include "GaudiKernel/ObjectFactory.h"
#include "GaudiKernel/ObjectVector.h"
#include "GaudiKernel/ObjectList.h"

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
// Object factory implementation for objects of class AcdDigi 
// ====================================================================
#include "GlastEvent/Digi/AcdDigi.h"
_ImplementHitContainedFactories(AcdDigi)

void DigiInstantiation()  {
    DLL_DECL_CONTAINEDOBJECTFACTORY( AcdDigi );
}
