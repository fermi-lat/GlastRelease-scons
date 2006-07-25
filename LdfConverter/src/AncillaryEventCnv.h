#ifndef AncillaryEventCnv_H
#define AncillaryEventCnv_H 1

#include "LdfBaseCnv.h"

//#include "Event/TopLevel/DigiEvent.h"
#include "AdfEvent/AdfEvent.h"

// Abstract factory to create the converter
template <class TYPE> class CnvFactory;


/** @class AncillaryEventCnv
 * @brief Concrete converter for the AncillaryEvent stored in the TDS 
 * /Event/AncillaryEvent
 *
 * $Header$
 */ 

class AncillaryEventCnv : public LdfBaseCnv { 

  friend class CnvFactory<AncillaryEventCnv>;

public: 
  static const CLID& classID()   
  {
    return CLID_AncillaryEvent;
  }

protected:

  AncillaryEventCnv(ISvcLocator* svc);

  virtual ~AncillaryEventCnv() { };

  /// override the EbfBaseCnv version to handle the conversion
  virtual StatusCode createObj(IOpaqueAddress* pAddress, DataObject*& refpObject);

  /// override the EbfBaseCnv version
  virtual StatusCode updateObj(int* data, DataObject* pObject);


};


#endif // AncillaryEventCnv_H
