#ifndef LdfDiagnosticCnv_H
#define LdfDiagnosticCnv_H 1

#include "LDfBaseCnv.h"
namespace LdfEvent { class DiagnosticData; }

extern const CLID& CLID_LdfDiagnosticData;

// Abstract factory to create the converter
template <class TYPE> class CnvFactory;

/** @class LdfDiagnosticCnv
 * @brief Concrete converter for the EM Diagnostic data stored in the TDS /Event/Diagnostic
 *
 * $Header$
 */ 

class LdfDiagnosticCnv : public LdfBaseCnv { 

  friend class CnvFactory<LdfDiagnosticCnv>;

public: 
  static const CLID& classID()   
  {
    return CLID_LdfDiagnosticData;
  }

protected:

  LdfDiagnosticCnv(ISvcLocator* svc);

  virtual ~LdfDiagnosticCnv() { };

  /// override the LdfBaseCnv version to handle the conversion
  virtual StatusCode createObj(IOpaqueAddress* pAddress, DataObject*& refpObject);

  /// override the LdfBaseCnv version
  virtual StatusCode updateObj(int* data, LdfEvent::DiagnosticData* pObject);


};


#endif // LdfDiagnosticCnv_H
