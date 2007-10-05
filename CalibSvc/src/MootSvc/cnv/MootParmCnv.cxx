// $Header$

/** @class MootParmCnv

  Converter from xml to TCDS CalAsym class

  @author J. Bogart
*/

#include "MootBaseCnv.h"
#include "mootCore/MootQuery.h"
#include "mootCore/FileDescrip.h"
#include "GaudiKernel/CnvFactory.h"
#include "CalibSvc/MootTds.h"
#include <string>
#include <vector>

// #include  something like the following:
// might not need them all
#include "GaudiKernel/IOpaqueAddress.h"
#include "GaudiKernel/DataObject.h"
#include "GaudiKernel/IAddressCreator.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/IConversionSvc.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/GenericAddress.h"


template <class TYPE> class CnvFactory;

class MootParmCnv : public MootBaseCnv {

  /// Friend needed for instantiation
  friend class CnvFactory<MootParmCnv>;
public:
  const CLID& objType() const;
  static const CLID& classID();
protected:

  MootParmCnv(ISvcLocator* svcs);

  virtual ~MootParmCnv() {}       // most likely nothing to do 

  // Reimplement
  virtual StatusCode createObj(IOpaqueAddress* addr, DataObject*& refpObject);


  virtual StatusCode updateObj(IOpaqueAddress* addr, DataObject*& refpObject);

private:

  /// Get info for latc source files assoc. with current LATC master key
  /// (Do nothing if it matches old cached version of key)
  /// This routine should return failure if subtype != MOOTSUBTYPE_latcParm
  //  StatusCode updateLatcParms(MootParmCol* parmCol);  TO BE WRITTEN

  StatusCode storeLatcParms(unsigned master, MootParmCol* parmCol);

  unsigned m_key;
};

static CnvFactory<MootParmCnv> s_factory;
const  ICnvFactory& MootParmCnvFactory = s_factory;

// Now start implementation

MootParmCnv::MootParmCnv(ISvcLocator* svc) :
  MootBaseCnv(svc, CLID_MootParm) {}

const CLID& MootParmCnv::objType() const {
  return CLID_MootParm;
}

const CLID& MootParmCnv::classID()  {
  return CLID_MootParm;
}

StatusCode MootParmCnv::createObj(IOpaqueAddress* addr,
                                  DataObject*& refpObject) {
  MOOTSUBTYPE sub = (MOOTSUBTYPE) addr->ipar()[1];
  unsigned newKey = addr->ipar()[2];

  MootParmCol* pParmCol = new MootParmCol(sub, 0);
  refpObject = pParmCol;
  switch (sub) {
  case MOOTSUBTYPE_latcParm:
    return storeLatcParms(newKey, pParmCol);
  default:   // shouldn't happen
    if (!m_log) m_log = new MsgStream(msgSvc(), "MootParmCnv");
    (*m_log) << MSG::ERROR << "bad subtype " << sub << endreq;
  }
  return StatusCode::SUCCESS;
}

StatusCode MootParmCnv::updateObj(IOpaqueAddress* addr,
                                  DataObject*& refpObject) {
  MootParmCol* pParmCol = dynamic_cast<MootParmCol*> (refpObject);
  unsigned newKey = addr->ipar()[2];
  if (newKey == pParmCol->m_key) return StatusCode::SUCCESS;

  MOOTSUBTYPE sub = (MOOTSUBTYPE) addr->ipar()[1];
  switch (sub) { // so far only one possibility
  case MOOTSUBTYPE_latcParm:
    return storeLatcParms(newKey, pParmCol);
  default:   // shouldn't happen
    if (!m_log) m_log = new MsgStream(msgSvc(), "MootParmCnv");
    (*m_log) << MSG::ERROR << "bad subtype " << sub << endreq;
  }
}

// Store parm information corresponding to key.  If successful, update
// local copy of key
StatusCode MootParmCnv::storeLatcParms(unsigned master, MootParmCol* parmCol) {
  if (!master) {
    // complain
    return StatusCode::FAILURE;
  }
  std::vector<MOOT::ParmOffline> poff;

  if (!m_q->getParmsFromMaster(master, poff)) return StatusCode::FAILURE;

  // Sort into order we need; store info
  //             ** TO-DO - BUNCH OF STUFF   **

}
