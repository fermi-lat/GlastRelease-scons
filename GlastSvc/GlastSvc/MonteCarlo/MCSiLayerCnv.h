// $Header$
#ifndef MCSiLayerCnv_H 
#define MCSiLayerCnv_H 1

// Include files
//#include "Gaudi/Kernel/Converter.h"
#include "src/EventSelector/BaseCnv.h"

// Abstract factory to create the converter
template <class TYPE> class CnvFactory;

//extern unsigned char SICB_StorageType;


//------------------------------------------------------------------------------
//
// ClassName:   MCSiLayerCnv
//  
// Description:  handles the conversion of TKR data stored as vectors of SiLayers
//
// Author:   
//
//------------------------------------------------------------------------------
class MCSiLayerCnv : public BaseCnv   {
  friend class CnvFactory<BaseCnv>;

public:
  /// Standard Constructor
    MCSiLayerCnv( ISvcLocator* svc );
  /// Standard Destructor
  virtual ~MCSiLayerCnv();
  //! access class ID
  static const CLID& classID();
  //! Update complete object
  virtual StatusCode updateObj(IOpaqueAddress* pAddress, DataObject* pObject);
  //! Create object
  virtual StatusCode createObj(IOpaqueAddress* pAddress, DataObject*& refpObject);

protected:

};
#endif // SiLayerCnv_H
