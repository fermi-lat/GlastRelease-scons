#ifndef IrfTkrLayerCnv_H 
#define IrfTkrLayerCnv_H 1

// Include files
//#include "Gaudi/Kernel/Converter.h"
#include "src/EventSelector/BaseCnv.h"

// Abstract factory to create the converter
template <class TYPE> class CnvFactory;


//------------------------------------------------------------------------------
//
// ClassName:   IrfTkrLayerCnv
//  
// Description:  handles the conversion of TKR data stored as vectors of TkrLayers
//
// Author:   
//
//------------------------------------------------------------------------------
class IrfTkrLayerCnv : public BaseCnv   {
  friend class CnvFactory<BaseCnv>;

public:
  /// Standard Constructor
    IrfTkrLayerCnv( ISvcLocator* svc );
  /// Standard Destructor
  virtual ~IrfTkrLayerCnv();
  //! access class ID
  static const CLID& classID();
  //! Update complete object
  virtual StatusCode updateObj(IOpaqueAddress* pAddress, DataObject* pObject);
  //! Create object
  virtual StatusCode createObj(IOpaqueAddress* pAddress, DataObject*& refpObject);

protected:

};
#endif // SiLayerCnv_H
