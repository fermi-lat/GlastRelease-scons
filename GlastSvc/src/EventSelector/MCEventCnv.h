// $Header$
#ifndef CNV_MCEVENTCNV_H 
#define CNV_MCEVENTCNV_H 1

// Include files
#include "ItemCnv.h"

// Forward declarations
class IDataDirectory;
class ISvcLocator;
class MCEvent;
// Abstract factory to create the converter
template <class T> class CnvFactory;

//------------------------------------------------------------------------------
//
// ClassName :    MCEventCnv
//  
// Description :  GLAST converter of Monte Carlo event header
//
// Author :       Markus Frank
//
//------------------------------------------------------------------------------
class MCEventCnv : public ItemCnv<MCEvent>      {
  friend class CnvFactory<MCEventCnv>;
protected:
  /// Standard Constructor
  MCEventCnv(ISvcLocator* svc);
  /// Standard Destructor
  virtual ~MCEventCnv();
  /// Update object from somewhere
  virtual StatusCode updateObj(int* data, MCEvent* pObject);
};
#endif // CNV_MCEVENTCNV_H
