#ifndef CNV_IrfEventCNV_H 
#define CNV_IrfEventCNV_H 1

// Include files
#include "ItemCnv.h"

// Forward declarations
class IDataDirectory;
class ISvcLocator;
class IrfEvent;
// Abstract factory to create the converter
template <class T> class CnvFactory;

//------------------------------------------------------------------------------
//
// ClassName :    IrfEventCnv
//  
// Description :  GLAST converter of Irf event header
//
//
//------------------------------------------------------------------------------
class IrfEventCnv : public ItemCnv<IrfEvent>      {
  friend class CnvFactory<IrfEventCnv>;
protected:
  /// Standard Constructor
  IrfEventCnv(ISvcLocator* svc);
  /// Standard Destructor
  virtual ~IrfEventCnv();
  /// Update object from somewhere
  virtual StatusCode updateObj(int* data, IrfEvent* pObject);
};
#endif // CNV_IrfEventCNV_H
