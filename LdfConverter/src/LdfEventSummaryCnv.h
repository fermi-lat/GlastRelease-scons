#ifndef LDFEVENTCNV_H
#define LDFEVENTCNV_H 1

#include "LdfBaseCnv.h"
namespace LdfEvent { class EventSummaryData; }

extern const CLID& CLID_LdfEventSummaryData;

// Abstract factory to create the converter
template <class TYPE> class CnvFactory;


/** @class LdfEventSummaryCnv
 * @brief Concrete converter for the Event header stored in the TDS /Event
 *
 * $Header$
 */ 

class LdfEventSummaryCnv : public LdfBaseCnv { 

  friend class CnvFactory<LdfEventSummaryCnv>;

public: 
  static const CLID& classID()   
  {
    return CLID_LdfEventSummaryData;
  }

protected:

  LdfEventSummaryCnv(ISvcLocator* svc);

  virtual ~LdfEventSummaryCnv() { };

  /// override the EbfBaseCnv version to handle the conversion
  virtual StatusCode createObj(IOpaqueAddress* pAddress, DataObject*& refpObject);

  /// override the EbfBaseCnv version
  virtual StatusCode updateObj(int* data, LdfEvent::EventSummaryData* pObject);


};


#endif // LdfEventCNV_H
