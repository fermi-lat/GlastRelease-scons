#ifndef IrfCalHitCnv_H 
#define IrfCalHitCnv_H 1

// Include files
//#include "Gaudi/Kernel/Converter.h"
#include "src/EventSelector/BaseCnv.h"

// Abstract factory to create the converter
template <class TYPE> class CnvFactory;

//------------------------------------------------------------------------------
//
// ClassName:   IrfCalHitCnv
//  
// Description: Gaudi Converter for IrfCalHitCnv
//
// Author:      
//
//------------------------------------------------------------------------------
class IrfCalHitCnv : public BaseCnv   {
  friend class CnvFactory<BaseCnv>;

public:
  /// Standard Constructor
  IrfCalHitCnv( ISvcLocator* svc );
  /// Standard Destructor
  virtual ~IrfCalHitCnv();
  /// class ID
  static const CLID& classID();
  /// Update complete object
  virtual StatusCode updateObj(IOpaqueAddress* pAddress, DataObject* pObject);
  /// Create object
  virtual StatusCode createObj(IOpaqueAddress* pAddress, DataObject*& refpObject);

protected:

};
#endif // IrfCalHitCnv_H
