// $Header$
#ifndef EVENTCNV_H
#define EVENTCNV_H 1


// Include files
#include "ItemCnv.h"


// Forward declarations
class Event;

// Abstract factory to create the converter
template <class TYPE> class CnvFactory;


//------------------------------------------------------------------------------
// (original stuff)
// ClassName :    EventCnv
//  
// Description :  converter of event header
//
// 
//
//------------------------------------------------------------------------------


class EventCnv : public ItemCnv<Event>	    {

  friend class CnvFactory<EventCnv>;

protected:

  /// Standard Constructor
  EventCnv(ISvcLocator* svc);

  /// Standard Destructor
  virtual ~EventCnv();

  virtual StatusCode updateObj(int* data, Event* pObject);

  //virtual StatusCode updateObjRefs(int* data, Event* pObject);

};


#endif // EVENTCNV_H
