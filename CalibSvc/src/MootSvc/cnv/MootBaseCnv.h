// $Header$
#ifndef CalibData_MootBaseCnv_h
#define CalibData_MootBaseCnv_h

/** @class MootBaseCnv 

  Base class for calibration converters from XML files to TCDS.
  All such converters need to do certain things, which are
  handled here.

  @author J. Bogart
*/
#include <string>
#include <vector>
#include "GaudiKernel/Converter.h"
#include "GaudiKernel/CnvFactory.h"
#include "GaudiKernel/MsgStream.h"

class ISvcLocator;
class GenericAddress;
class IMootSvc;

namespace MOOT {
  class MootQuery;
}

class MootBaseCnv : public Converter {

public:
  /**
     Constructor for this converter
     @param svc a ISvcLocator interface to find services
     @param clid the type of object the converter is able to convert
   */
  MootBaseCnv(ISvcLocator* svc, const CLID& clid);

  // Needed to satisay IConverter interface
  /// Retrieve the class type of the data store the converter uses.
  virtual long repSvcType() const {return Converter::i_repSvcType();}
  
  virtual ~MootBaseCnv() {};

  virtual StatusCode initialize();

  virtual StatusCode finalize();

  /**
   Create the transient representation of an object, given an opaque
   address.  All MootBaseCnv does is set the type and subtype enums
  */
  virtual StatusCode createObj(IOpaqueAddress* addr,
                               DataObject*& refpObject);

  // Nothing for MootBaseCnv to do at update time, so don't bother to
  // reimplement

  IMootSvc* getMootSvc() {
    return m_mootSvc;
  }


  static const unsigned char storageType();

protected:

  IMootSvc* m_mootSvc;
  MOOT::MootQuery* m_q;
  MsgStream*           m_log;

};
#endif
