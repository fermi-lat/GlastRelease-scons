//====================================================================
// IRFDbConverters.cpp
//--------------------------------------------------------------------
/*!
Defines the macros for the Root converters for the various entities contained
in the TDS.
Based upon MonteCarloConverters.cpp created by Markus Frank within DbCnv.
*/   

#define DBCNV_IRFDBCONVERTERS_CPP  1

// Include files
#include "GlastEvent/TopLevel/ObjectVector.h"
#include "GlastEvent/TopLevel/ObjectList.h"
#include "DbCnv/DbContainerConverter.h"


//#include "DbCnv/DbUserConverter.h"


// ====================================================================
// Converter implementation for objects of class MCCalorimeterHit
// ====================================================================
#include "GlastEvent/MonteCarlo/MCCalorimeterHit.h"
_ImplementContainerConverters(MCCalorimeterHit)

// ====================================================================
// Converter implementation for objects of class MCEcalFacePlaneHit
// ====================================================================
#include "GlastEvent/MonteCarlo/MCACDHit.h"
_ImplementContainerConverters(MCACDHit)
// alternative to this method
//static const DbUserCnvFactory < ObjectVector<MCACDHit> > s_MCACDHitConverterFactory;
//const ICnvFactory& DbMCACDHitConverterFactory = s_MCACDHitConverterFactory;

// ====================================================================
// Converter implementation for objects of class MCRichPhotoDetectorHit
// ====================================================================
#include "GlastEvent/MonteCarlo/MCTKRHit.h"
_ImplementContainerConverters(MCTKRHit)


// ====================================================================
// Converter implementation for objects of class MCHitBase
// Objects cannot be instantiated!
// ===================================================================
//_ImplementContainerConverters(MCHitBase)
