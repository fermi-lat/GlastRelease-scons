
// $Header$
#ifndef MCCalorimeterHitCnv_H 
#define MCCalorimeterHitCnv_H 1

// Include files
//#include "Gaudi/Kernel/Converter.h"
#include "src/EventSelector/BaseCnv.h"

// Abstract factory to create the converter
template <class TYPE> class CnvFactory;

//------------------------------------------------------------------------------
//
// ClassName:   MCCalorimeterHitCnv
//  
// Description: Gaudi Converter for MCCalorimeterHitCnv
//
// Author:      
//
//------------------------------------------------------------------------------
class MCCalorimeterHitCnv : public BaseCnv   {
  friend class CnvFactory<BaseCnv>;

public:
  /// Standard Constructor
  MCCalorimeterHitCnv( ISvcLocator* svc );
  /// Standard Destructor
  virtual ~MCCalorimeterHitCnv();
  /// class ID
  static const CLID& classID();
  /// Update complete object
  virtual StatusCode updateObj(IOpaqueAddress* pAddress, DataObject* pObject);
  /// Create object
  virtual StatusCode createObj(IOpaqueAddress* pAddress, DataObject*& refpObject);

protected:

};
#endif // MCCalorimeterHitCnv_H
