/**
  @file MetadataCnv.h

  @author Joanne Bogart
  $Header$
*/

#ifndef MetdataCnv_H
#define MetdataCnv_H

#include "GaudiKernel/Converter.h"

class ICalibCnvSvc;

template <class TYPE> class CnvFactory;

namespace calibUtil {
  class Metadata;
}

/** @class MetadataCnv
  Class definition for trivial calib. data converter: reads in
  metadata.
*/
class MetadataCnv : public Converter {
  // needs to call our protected constructor
  friend class CnvFactory<MetadataCnv>;  

public:
  virtual StatusCode initialize();

  virtual StatusCode finalize();

  /**
   * Creates the transient representation of an object.
   * @param addr the address of the object representation
   * @param refpObject the object created
   * @return status depending on the completion of the call
   */
  virtual StatusCode createObj (IOpaqueAddress* addr,
                                DataObject*& refpObject);
  /**
   * Updates the transient object from the other representation.
   * @param pAddress the address of the object representation
   * @param pObject the object updated
   *  @return status depending on the completion of the call
   */
  virtual StatusCode updateObj (IOpaqueAddress* pAddress,
                                DataObject* pObject);
  
  // Don't override createRep or updateRep.  (Base class just returns
  // success.)  For now we don't support TDDS --> persistent 

protected:

  /**
     Constructor is protected.  The converter can only be constructed
     by its friend, the factory.
  */

  MetadataCnv(ISvcLocator* svc);

  virtual ~MetadataCnv() {}

  ICalibCnvSvc* m_CalibCnvSvc;
  IDetDataSvc* m_DetDataSvc;
}

#endif
