//  ====================================================================
//  BaseCnv.h
//  --------------------------------------------------------------------
//
//  Package   : GlastSvc/src/EventSelector
//
//
//  ====================================================================
#ifndef _BASECNV_H
#define _BASECNV_H 1

// Include files
#include <vector>
#include <string>
#include "Gaudi/Kernel/Converter.h"
#include "IEventCnvSvc.h"

// HMA:  Took SicbBaseCnv and just renamed the class...

// forward declarations
class IEventCnvSvc;
class IGlastDetSvc;

// external declarations
extern unsigned char SICB_StorageType;

/**

    Description:

*/
class BaseCnv : public Converter      {
protected:
  /// Pointer to (enhanced) event conversion service interface
  IEventCnvSvc*                     m_CnvSvc;
  /// Leaf objects giving the locations of the objects in the data store
  std::vector<IEventCnvSvc::Leaf>   m_leaves;
  // access to the GlastDetSvc
  IGlastDetSvc*       m_detSvc;


public:
  /// Standard Constructor
  BaseCnv(const CLID& clid, ISvcLocator* svc);

  /// Standard Destructor
  virtual ~BaseCnv();

  /// Storage type and class ID
  static const unsigned char storageType()      { 
    return SICB_StorageType; 
  }
  /// Initialization
  virtual StatusCode initialize();

  /// finalization
  virtual StatusCode finalize();

  /// Declare target leaf
  void declareObject(const std::string& fullPath, const CLID& clid, const std::string& bank, long par=0);

  /// Convert the transient object to the requested representation.
  virtual StatusCode createRep(DataObject* pObject, IOpaqueAddress*& refpAddress);

  /// Resolve the references of the converted object. 
  virtual StatusCode fillRepRefs(IOpaqueAddress* pAddress,DataObject* pObject);

  /// Update the converted representation of a transient object.
  virtual StatusCode updateRep(IOpaqueAddress* pAddress, DataObject* pObject);

  /// Update the references of an already converted object.
  virtual StatusCode updateRepRefs(IOpaqueAddress* pAddress, DataObject* pObject);
};


#endif // _BASECNV_H
