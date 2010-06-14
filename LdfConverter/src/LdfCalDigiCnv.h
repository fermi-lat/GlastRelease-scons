#ifndef LdfCalDigiCnv_H
#define LdfCalDigiCnv_H 1

#include "LdfBaseCnv.h"

namespace Event{ class EventHeader; }
extern const CLID& CLID_CalDigi;

// Abstract factory to create the converter
template <class TYPE> class CnvFactory;


/** @class LdfCalDigiCnv
 * @brief Concrete converter for the Event header stored in the TDS /Event
 *
 * $Header$
 */ 

class LdfCalDigiCnv : public LdfBaseCnv { 

  friend class CnvFactory<LdfCalDigiCnv>;

public: 
  static const CLID& classID()   
  {
    return CLID_CalDigi;
  }

protected:

  LdfCalDigiCnv(ISvcLocator* svc);

  virtual ~LdfCalDigiCnv() { };

  /// override the EbfBaseCnv version to handle the conversion
  virtual StatusCode createObj(IOpaqueAddress* pAddress, DataObject*& refpObject);

  /// override the LdfBaseCnv version
  virtual StatusCode updateObj(int* data, Event::EventHeader* pObject);


};


#endif // LdfCalDigiCnv_H
