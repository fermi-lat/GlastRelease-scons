// LdGlastDataCnv.h: interface for the LdGlastDataCnv class.
//
//////////////////////////////////////////////////////////////////////

#ifndef LdGlastDataCnv_H 
#define LdGlastDataCnv_H 1

// Include files
//#include "Gaudi/Kernel/Converter.h"
#include "src/EventSelector/BaseCnv.h"

extern const CLID& CLID_LdGlastDataCnv;

// Abstract factory to create the converter
template <class TYPE> class CnvFactory;



// Original Author: Ian Gable (03/02/2001)

/*! Handles the conversion of LdGlastData class currently uses the IRF
    information.
*/

class LdGlastDataCnv : public BaseCnv   {
  friend class CnvFactory<BaseCnv>;

public:
  /// Standard Constructor
  LdGlastDataCnv( ISvcLocator* svc );
  /// Standard Destructor
  virtual ~LdGlastDataCnv();
  /// class ID
  static const CLID& classID();
  /// Update complete object
  virtual StatusCode updateObj(IOpaqueAddress* pAddress, DataObject* pObject);
  /// Create object
  virtual StatusCode createObj(IOpaqueAddress* pAddress, DataObject*& refpObject);

protected:

};
#endif // LdGlastDataCnv_H
