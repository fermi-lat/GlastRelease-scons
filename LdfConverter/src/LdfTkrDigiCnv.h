#ifndef LdfTkrDigiCnv_H
#define LdfTkrDigiCnv_H 1

#include "LdfBaseCnv.h"

#include "Event/Digi/TkrDigi.h"

// Abstract factory to create the converter
template <class TYPE> class CnvFactory;


/** @class LdfTkrDigiCnv
 * @brief Concrete converter for the Event header stored in the TDS /Event
 *
 * $Header$
 */ 

class LdfTkrDigiCnv : public LdfBaseCnv { 

  friend class CnvFactory<LdfTkrDigiCnv>;

public: 
  static const CLID& classID()   
  {
    return CLID_TkrDigi;
  }

protected:

  LdfTkrDigiCnv(ISvcLocator* svc);

  virtual ~LdfTkrDigiCnv() { };

  /// override the EbfBaseCnv version to handle the conversion
  virtual StatusCode createObj(IOpaqueAddress* pAddress, DataObject*& refpObject);

  /// override the EbfBaseCnv version
  virtual StatusCode updateObj(int* data, Event::TkrDigi* pObject);


};


#endif // LdfTkrDigiCnv_H
