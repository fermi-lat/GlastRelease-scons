#ifndef EVENTCNV_H
#define EVENTCNV_H 1

#include "BaseCnv.h"

class Event;
extern const CLID& CLID_Event;

// Abstract factory to create the converter
template <class TYPE> class CnvFactory;


/** @class EventCnv
 * @brief Concrete Converter for the Event header
 *
 * $Header$
 */ 

class EventCnv : public BaseCnv { //public ItemCnv<Event>	    {

  friend class CnvFactory<EventCnv>;

public: 
  static const CLID& classID()   
  {
    return CLID_Event;
  }

protected:

  EventCnv(ISvcLocator* svc);

  virtual ~EventCnv() { };

  virtual StatusCode createObj(IOpaqueAddress* pAddress, DataObject*& refpObject);

  virtual StatusCode updateObj(int* data, Event* pObject);


};


#endif // EVENTCNV_H
