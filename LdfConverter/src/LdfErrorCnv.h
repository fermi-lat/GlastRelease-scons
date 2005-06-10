#ifndef LdfErrorCnv_H
#define LdfErrorCnv_H 1

#include "LdfBaseCnv.h"
#include "LdfEvent/ErrorData.h"

// Abstract factory to create the converter
template <class TYPE> class CnvFactory;

/** @class LdfErrorCnv
 * @brief Concrete converter for the error data stored in the TDS /Event/Diagnostic
 *
 * $Header$
 */ 

class LdfErrorCnv : public LdfBaseCnv { 

  friend class CnvFactory<LdfErrorCnv>;

public: 
  static const CLID& classID()   
  {
    return CLID_LdfErrorData;
  }

protected:

  LdfErrorCnv(ISvcLocator* svc);

  virtual ~LdfErrorCnv() { };

  /// override the LdfBaseCnv version to handle the conversion
  virtual StatusCode createObj(IOpaqueAddress* pAddress, DataObject*& refpObject);

  /// override the LdfBaseCnv version
  virtual StatusCode updateObj(int* data, LdfEvent::ErrorData* pObject);


};


#endif // LdfErrorCnv_H
