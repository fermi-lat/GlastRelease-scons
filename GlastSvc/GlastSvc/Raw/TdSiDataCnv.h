// TdSiDataCnv.h: interface for the TdSiDataCnv class.
//
//////////////////////////////////////////////////////////////////////

#ifndef TdSiDataCnv_H 
#define TdSiDataCnv_H 1

// Include files
//#include "Gaudi/Kernel/Converter.h"
#include "src/EventSelector/BaseCnv.h"

extern const CLID& CLID_TdSiDataCnv;

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
class TdSiDataCnv : public BaseCnv   {
  friend class CnvFactory<BaseCnv>;

public:
  /// Standard Constructor
  TdSiDataCnv( ISvcLocator* svc );
  /// Standard Destructor
  virtual ~TdSiDataCnv();
  /// class ID
  static const CLID& classID();
  /// Update complete object
  virtual StatusCode updateObj(IOpaqueAddress* pAddress, DataObject* pObject);
  /// Create object
  virtual StatusCode createObj(IOpaqueAddress* pAddress, DataObject*& refpObject);

protected:

};
#endif // TdSiDataCnv_H
