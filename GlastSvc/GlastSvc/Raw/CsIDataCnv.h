// CsIDataCnv.h: interface for the CsIDataCnv class.
//
//////////////////////////////////////////////////////////////////////

#ifndef CsIDataCnv_H 
#define CsIDataCnv_H 1

// Include files
//#include "Gaudi/Kernel/Converter.h"
#include "src/EventSelector/BaseCnv.h"

extern const CLID& CLID_CsIDataCnv;

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
class CsIDataCnv : public BaseCnv   {
  friend class CnvFactory<BaseCnv>;

public:
  /// Standard Constructor
  CsIDataCnv( ISvcLocator* svc );
  /// Standard Destructor
  virtual ~CsIDataCnv();
  /// class ID
  static const CLID& classID();
  /// Update complete object
  virtual StatusCode updateObj(IOpaqueAddress* pAddress, DataObject* pObject);
  /// Create object
  virtual StatusCode createObj(IOpaqueAddress* pAddress, DataObject*& refpObject);

protected:

};
#endif // CsIDataCnv_H
