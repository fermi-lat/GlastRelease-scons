#ifndef LdfDigiCnv_H
#define LdfDigiCnv_H 1

#include "LdfBaseCnv.h"

namespace Event{ class EventHeader; }
extern const CLID& CLID_DigiEvent;

// Abstract factory to create the converter
template <class TYPE> class CnvFactory;


/** @class LdfDigiCnv
 * @brief Concrete converter for the Event header stored in the TDS /Event
 *
 * $Header$
 */ 

class LdfDigiCnv : public LdfBaseCnv { 

  friend class CnvFactory<LdfDigiCnv>;

public: 
  static const CLID& classID()   
  {
    return CLID_DigiEvent;
  }

protected:

  LdfDigiCnv(ISvcLocator* svc);

  virtual ~LdfDigiCnv() { };

  /// override the EbfBaseCnv version to handle the conversion
  virtual StatusCode createObj(IOpaqueAddress* pAddress, DataObject*& refpObject);

  /// override the EbfBaseCnv version
  virtual StatusCode updateObj(int* data, Event::EventHeader* pObject);


};


#endif // LdfDigiCnv_H
