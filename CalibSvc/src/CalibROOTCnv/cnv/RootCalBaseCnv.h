// $Header$
#ifndef RootCalBaseCnv_h
#define RootCalBaseCnv_h

/** @class RootCalBaseCnv 

  Base class for CAL calibration converters from ROOT files to TCDS.
  All such converters need to do certain things, which are
  handled here.  Methods common to *all* calibrations are in the
  base class RootBaseCnv

  @author J. Bogart
*/

#include "RootBaseCnv.h"

namespace CalibData {
  class CalCalibBase;
}

class RootCalBaseCnv : public RootBaseCnv {

public:

  RootCalBaseCnv(ISvcLocator* svc, const CLID& clid);

  virtual ~RootCalBaseCnv() {};


  //  virtual StatusCode finalize();

  /**
   Create the transient representation of an object, given an opaque
   address.  This and the following update method comprise the core 
   functionality of calibration converters.
  */
  //  virtual StatusCode createObj(IOpaqueAddress* addr,
  //                               DataObject*& refpObject);

  //  virtual StatusCode i_createObj (const DOM_Element& element,
  //                                  DataObject*& refpObject);

  //  virtual StatusCode i_processObj(DataObject* pObject,
  //

  /// Convenience routine used by most CAL calibration types, which
  /// have a <dimension> element describing how the remainder of the
  /// data is laid out.  Read from TDS; store information internally 
  /// in protected members.  
  StatusCode readDimension(CalibData::CalCalibBase*  pCalBase);

  // Might need another one reading from Root class into protected members

  /*
    Not sure yet what the analogous thing to findFirstRange, findNextRange
    ought to be..likely don't need it at all since ROOT persistent form 
    comes fully labeled with CalXtalId
  */

protected:
  /**
     Given a pointer to a TDS object which can be cast to "our" type, fill
     in corresponding information in the corresponding root class

     @param pTDSObj   Pointer to tds object to be converted
     @param pRootObj  Pointer to destination root object

  */
  virtual StatusCode fillRoot(CalibData::CalibBase* pTDSObj, 
                              TObject* pRootObj);

  /**
     Read in object from specified branch. Don't need tree name; it's
     always Calib
  */
  virtual StatusCode readRootObj(const std::string& branch, TObject*& pCalib);

  /// Keep dimension info here
  unsigned m_nRow;
  unsigned m_nCol;
  unsigned m_nLayer;
  unsigned m_nXtal;
  unsigned m_nFace;
  unsigned m_nRange;
  unsigned m_nDacCol;  // not sure this is needed

};

#endif
