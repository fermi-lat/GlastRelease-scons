#ifndef IrfAcdHitCnv_H 
#define IrfAcdHitCnv_H 1

// Include files
//#include "Gaudi/Kernel/Converter.h"
#include "src/EventSelector/BaseCnv.h"

// Abstract factory to create the converter
template <class TYPE> class CnvFactory;


//------------------------------------------------------------------------------
//ClassName:   IrfAcdHitCnv
/*!
A converter class that handles ACD data from the IRF file
*/
// Author:   
//
//------------------------------------------------------------------------------
class IrfAcdHitCnv : public BaseCnv   {
  friend class CnvFactory<BaseCnv>;

public:
  //! Standard Constructor
  IrfAcdHitCnv( ISvcLocator* svc );
  //! Standard Destructor
  virtual ~IrfAcdHitCnv();
  //! access class ID
  static const CLID& classID();
  //! Update complete object
  virtual StatusCode updateObj(IOpaqueAddress* pAddress, DataObject* pObject);
  //! Create object
  virtual StatusCode createObj(IOpaqueAddress* pAddress, DataObject*& refpObject);

protected:

};
#endif // ACDhitCnv_H
