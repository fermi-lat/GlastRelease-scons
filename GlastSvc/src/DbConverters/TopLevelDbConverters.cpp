//	====================================================================
//	TopLevelDbConverters.cpp
//	--------------------------------------------------------------------
//
/*!
Defines the macros for the Root converters for the various entities contained
in the TDS.
Based upon TopLevelConverters.cpp created by Markus Frank within DbCnv.
*/   

#define DBCNV_TOPLEVELCONVERTERS_CPP  1

// Include files
#include "DbCnv/DbGenericConverter.h"
#include "GlastEvent/TopLevel/ObjectVector.h"
#include "GlastEvent/TopLevel/ObjectList.h"
#include "DbCnv/DbContainerConverter.h"

// ====================================================================
// Converter implementation for objects of class DataObject
// ===================================================================
#include "Gaudi/Kernel/DataObject.h"
//_ImplementDataObjectCreator(DataObject)
_ImplementConverter(DataObject)

// ====================================================================
// Converter implementation for objects of class Run
// ===================================================================
//#include "GlastEvent/TopLevel/Run.h"
//_ImplementConverter(Run)

// ====================================================================
// Converter implementation for objects of class Event
// ===================================================================
#include "GlastEvent/TopLevel/Event.h"
_ImplementConverter(Event)

// ====================================================================
// Converter implementation for objects of class MCEvent
// ===================================================================
#include "GlastEvent/TopLevel/MCEvent.h"
_ImplementConverter(MCEvent)

// ====================================================================
// Converter implementation for objects of class IrfEvent
// ===================================================================
#include "GlastEvent/TopLevel/IrfEvent.h"
_ImplementConverter(IrfEvent)

// ====================================================================
// Converter implementation for objects of class ContainedObject
// ====================================================================
#include "DbCnv/DbContainerConverter.h"
#include "Gaudi/Kernel/ContainedObject.h"
_ImplementContainerConverters(ContainedObject)
