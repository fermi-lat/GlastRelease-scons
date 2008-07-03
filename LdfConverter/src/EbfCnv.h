#ifndef EbfCnv_H
#define EbfCnv_H 1

#include "LdfBaseCnv.h"
#include "EbfWriter/Ebf.h"

// Abstract factory to create the converter
template <class TYPE> class CnvFactory;


/** @class EbfCnv
 * @brief Concrete converter for the EBF event 
 *
 * $Header$
 */ 

class EbfCnv : public LdfBaseCnv { 

  friend class CnvFactory<EbfCnv>;

public: 
  static const CLID& classID()   
  {
    return CLID_Ebf;
  }

protected:

  EbfCnv(ISvcLocator* svc);

  virtual ~EbfCnv() { };

  /// override the LdfBaseCnv version to handle the conversion
  virtual StatusCode createObj(IOpaqueAddress* pAddress, DataObject*& refpObject);

  /// override the LdfBaseCnv version
  virtual StatusCode updateObj(int* data, EbfWriterTds::Ebf* pObject);



};


#endif // EbfCnv_H
