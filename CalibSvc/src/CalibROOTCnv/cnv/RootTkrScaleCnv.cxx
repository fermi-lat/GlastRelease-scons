// $Header$

#include <string>
#include <ios>
#include "TTree.h"
#include "TFile.h"
#include "TBranch.h"

#include "GaudiKernel/CnvFactory.h"
#include "GaudiKernel/IOpaqueAddress.h"
#include "GaudiKernel/DataObject.h"
#include "GaudiKernel/IAddressCreator.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/IConversionSvc.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/GenericAddress.h"

#include "CalibSvc/ICalibRootSvc.h"      // maybe
#include "CalibSvc/ICalibMetaCnvSvc.h"

#include "RootTkrScaleCnv.h"

#include "CalibData/Tkr/TkrScale.h"
// #include "CalibData/CalibTime.h"
#include "commonRootData/idents/TkrId.h"
#include "calibRootData/Tkr/ChargeScale.h"
#include "idents/TkrId.h"

// Temporary.  Hope to find a better way to do this
#include "CalibData/CalibModel.h"

//static CnvFactory<RootTkrScaleCnv> s_factory;
//const  ICnvFactory& RootTkrScaleCnvFactory = s_factory;
DECLARE_CONVERTER_FACTORY(RootTkrScaleCnv);

RootTkrScaleCnv::RootTkrScaleCnv( ISvcLocator* svc) :
  RootTkrBaseCnv(svc, CLID_Calib_TKR_ChargeScale) { 
}


const CLID& RootTkrScaleCnv::objType() const {
  return CLID_Calib_TKR_ChargeScale;
}

const CLID& RootTkrScaleCnv::classID() {
  return CLID_Calib_TKR_ChargeScale;
}

StatusCode RootTkrScaleCnv::i_createObj (const std::string& fname,
                                       DataObject*& refpObject) {
  using CalibData::TkrBase;
  using CalibData::TkrScaleCol;
  // open file
  openRead(fname);

  // Call CalibData::TkrScaleCol constructor
  // The next two executable lines are the only ones in this routine which
  // are unique to this calibration type.  Remaining code specific to
  // this calib type is in readUnis
  TkrScaleCol* scaleCol = new TkrScaleCol();
  refpObject = scaleCol;

  setBaseInfo(scaleCol);
  StatusCode ret;
  for (unsigned iTow = 0; iTow < TKRBASE_MAXTOWER; iTow++) {
    TTree* tree = findTower(iTow);
    if (tree) {
      // handle generic tracker part
      ret = readTower(tree, iTow, scaleCol);
      if (ret != StatusCode::SUCCESS) {
        closeRead();
        return ret;
      }

      // read in Scale info for each uniplane
      ret = readUnis(tree, iTow,  scaleCol);
      if (ret != StatusCode::SUCCESS) {
        closeRead();
        return ret;
      }
    }
  }
  closeRead();
  return ret;
}

StatusCode RootTkrScaleCnv::readUnis(TTree* tree, int iTow,
                                   CalibData::TkrScaleCol* col) {
  using CalibData::TkrScaleCol;
  using CalibData::UniBase;
  using CalibData::TkrScaleObj;

  MsgStream log(msgSvc(), "RootTkrScaleCnv");

  // set branch 
  TBranch* branch = tree->GetBranch("calibRootData::ChargeScaleUnilayer");

  // find #unis
  Stat_t nEntries = branch->GetEntries();

  for (unsigned ix = 0; ix < nEntries; ix++) {  // process a unilayer

    calibRootData::ChargeScaleUnilayer* rootUni = 
      new calibRootData::ChargeScaleUnilayer();
    TObject* pObj = rootUni;

    StatusCode ret = 
      readRootObj(tree, "calibRootData::ChargeScaleUnilayer", pObj, ix);
    if (ret != StatusCode::SUCCESS) {
      log << MSG::ERROR << "Failed to read ChargeScale unilayer object" 
          << endreq;
      return ret;
    }
    idents::TkrId id;
  
    const commonRootData::TkrId& rootId = rootUni->getId();

    ret = convertId(rootId, &id);
    if (ret != StatusCode::SUCCESS) {
      log << MSG::ERROR << "Illegal TkrId in ChargeScale root calib. object" 
          << endreq;
      return ret;
    }

    if (!checkTower(id, iTow)) {
      // complain == uni doesn't belong to the right tower
      log << MSG::ERROR 
          << "Charge scale unilayer TkrId inconsistent with tower" << endreq;
      return StatusCode::FAILURE; 
    }

    int nObjs = rootUni->getNChildren();
    unsigned uniIx = 2*id.getTray() + id.getBotTop();

    std::vector<UniBase*> *unis = col->getUnis(iTow);
    CalibData::TkrScaleUni* dest;
    if (!(*unis)[uniIx]) {
      // have to make one
      (*unis)[uniIx]= makeUni(col);
    }
  
    dest = dynamic_cast<CalibData::TkrScaleUni*>((*unis)[uniIx]);
    if (!dest) {
      // tilt!
      log << MSG::ERROR << "Unable to create CalibData::TkrScaleUni object" 
          << endreq;
      return StatusCode::FAILURE;
    }  
    // For now, check if nObjs match.  If not, delete old 
    // (if exists) and allocate new.
    // But this isn't really quite right.  Can find out how
    // many strip objects are in ROOT object, but really what
    // we need to know is max strip id occurring.
    // Could add a "maxId"  field to calibRootData::ChargeScaleUnilayer
    dest->resize(nObjs, rootUni->childIsStrip() );

    if (nObjs > 0) {   // copy strip or gtfe info
      for (unsigned iObj = 0;  iObj < (unsigned) nObjs; iObj++) {
        const calibRootData::ChargeScaleObj* rObj = rootUni->getObj(iObj);
        TkrScaleObj tdsObj(rObj->getId(), rObj->getScale(),
                             rObj->getError(), rObj->getChi2(), rObj->getDf());
        // bool ok = 
        dest->putObj(tdsObj);
      }
    }    
    delete rootUni;  // is this required or forbidden?
    rootUni = 0;
  }    // end unilayer
  return StatusCode::SUCCESS;
}
