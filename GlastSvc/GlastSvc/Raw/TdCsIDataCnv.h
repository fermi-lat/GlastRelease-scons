// TdCsIDataCnv.h: interface for the TdCsIDataCnv class.
//
//////////////////////////////////////////////////////////////////////

#ifndef TdCsIDataCnv_H 
#define TdCsIDataCnv_H 1

// Include files
//#include "Gaudi/Kernel/Converter.h"
#include "src/EventSelector/BaseCnv.h"

extern const CLID& CLID_TdCsIDataCnv;

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
class TdCsIDataCnv : public BaseCnv   {
  friend class CnvFactory<BaseCnv>;

public:
  /// Standard Constructor
  TdCsIDataCnv( ISvcLocator* svc );
  /// Standard Destructor
  virtual ~TdCsIDataCnv();
  /// class ID
  static const CLID& classID();
  /// Update complete object
  virtual StatusCode updateObj(IOpaqueAddress* pAddress, DataObject* pObject);
  /// Create object
  virtual StatusCode createObj(IOpaqueAddress* pAddress, DataObject*& refpObject);

protected:

};
#endif // TdCsIDataCnv_H
