#ifndef EVENTCNV_H
#define EVENTCNV_H 1

#include "BaseCnv.h"

namespace Event{ class EventHeader; }
extern const CLID& CLID_Event;

// Abstract factory to create the converter
template <class TYPE> class CnvFactory;


/** @class EventCnv
 * @brief Concrete converter for the Event header stored in the TDS /Event
 *
 * $Header$
 */ 

class EventCnv : public BaseCnv { 

  friend class CnvFactory<EventCnv>;

public: 
  static const CLID& classID()   
  {
    return CLID_Event;
  }

protected:

  EventCnv(ISvcLocator* svc);

  virtual ~EventCnv() { };

  /// override the BaseCnv version to handle the conversion
  virtual StatusCode createObj(IOpaqueAddress* pAddress, DataObject*& refpObject);

  /// override the BaseCnv version
  virtual StatusCode updateObj(int* data, Event::EventHeader* pObject);


};


#endif // EVENTCNV_H
