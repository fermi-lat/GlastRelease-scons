// $Header$
#ifndef CNV_MCEVENTCNV_H 
#define CNV_MCEVENTCNV_H 1

// Include files
#include "BaseCnv.h"

// Forward declarations
class IDataDirectory;
class ISvcLocator;
class MCEvent;


// Abstract factory to create the converter
template <class T> class CnvFactory;

//------------------------------------------------------------------------------
//
// ClassName :    MCEventCnv
//  
// Description :  GLAST converter of Monte Carlo event header
//
// Author :      
//
//------------------------------------------------------------------------------
class MCEventCnv : public BaseCnv { // ItemCnv<MCEvent>      {
  friend class CnvFactory<MCEventCnv>;
public:
  /// Standard Constructor
  MCEventCnv(ISvcLocator* svc);
  /// Standard Destructor
  virtual ~MCEventCnv();
  static const CLID& classID();

  /// Create object
  virtual StatusCode createObj(IOpaqueAddress* pAddress, DataObject*& refpObject);

private:
};
#endif // CNV_MCEVENTCNV_H
