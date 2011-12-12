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

#include "RootTkrTotCnv.h"

#include "CalibData/Tkr/TkrTot.h"
// #include "CalibData/CalibTime.h"
#include "commonRootData/idents/TkrId.h"
#include "calibRootData/Tkr/Tot.h"
#include "idents/TkrId.h"

// Temporary.  Hope to find a better way to do this
#include "CalibData/CalibModel.h"

//static CnvFactory<RootTkrTotCnv> s_factory;
//const  ICnvFactory& RootTkrTotCnvFactory = s_factory;
DECLARE_CONVERTER_FACTORY(RootTkrTotCnv);

RootTkrTotCnv::RootTkrTotCnv( ISvcLocator* svc) :
  RootTkrBaseCnv(svc, CLID_Calib_TKR_TOTSignal) { 
}


const CLID& RootTkrTotCnv::objType() const {
  return CLID_Calib_TKR_TOTSignal;
}

const CLID& RootTkrTotCnv::classID() {
  return CLID_Calib_TKR_TOTSignal;
}

StatusCode RootTkrTotCnv::i_createObj (const std::string& fname,
                                       DataObject*& refpObject) {
  using CalibData::TkrBase;
  using CalibData::TkrTotCol;
  // open file
  openRead(fname);

  // Call CalibData::TkrTotCol constructor
  // The next two executable lines are the only ones in this routine which
  // are unique to this calibration type.  Remaining code specific to
  // this calib type is in readUnis
  TkrTotCol* totCol = new TkrTotCol();
  refpObject = totCol;

  setBaseInfo(totCol);
  StatusCode ret;
  for (unsigned iTow = 0; iTow < TKRBASE_MAXTOWER; iTow++) {
    TTree* tree = findTower(iTow);
    if (tree) {
      // handle generic tracker part
      ret = readTower(tree, iTow, totCol);
      if (ret != StatusCode::SUCCESS) {
        closeRead();
        return ret;
      }
      // read in Tot info for each uniplane
      ret = readUnis(tree, iTow,  totCol);
      if (ret != StatusCode::SUCCESS) {
        closeRead();
        return ret;
      }
    }
  }
  closeRead();
  return ret;
}

StatusCode RootTkrTotCnv::readUnis(TTree* tree, int iTow,
                                   CalibData::TkrTotCol* col) {
  using CalibData::TkrTotCol;
  using CalibData::UniBase;
  using CalibData::TkrTotStrip;

  MsgStream log(msgSvc(), "RootTkrTotCnv");

  // set branch 
  TBranch* branch = tree->GetBranch("calibRootData::TotUnilayer");

  // find #unis
  Stat_t nEntries = branch->GetEntries();

  for (unsigned ix = 0; ix < nEntries; ix++) {  // process a unilayer

    calibRootData::TotUnilayer* rootUni = new calibRootData::TotUnilayer();
    TObject* pObj = rootUni;

    StatusCode ret = 
      readRootObj(tree, "calibRootData::TotUnilayer", pObj, ix);
    if (ret != StatusCode::SUCCESS) {
      log << MSG::ERROR << "Failed to read ToT unilayer object" << endreq;
      return ret;
    }
    idents::TkrId id;
  
    const commonRootData::TkrId& rootId = rootUni->getId();

    ret = convertId(rootId, &id);
    if (ret != StatusCode::SUCCESS) {
      log << MSG::ERROR << "Illegal TkrId in ToT root calib. object" << endreq;
      return ret;
    }

    if (!checkTower(id, iTow)) {
      // complain == uni doesn't belong to the right tower
      log << MSG::ERROR << "ToT unilayer TkrId inconsistent with tower" 
          << endreq;
      return StatusCode::FAILURE; 
    }

    int nStrips = rootUni->getNStrips();
    unsigned uniIx = 2*id.getTray() + id.getBotTop();
    /*
    TkrTotUni* dest = dynamic_cast<TkrTotUni* > (&tdsTow->m_unis[uniIx]);
    */
    std::vector<UniBase*> *unis = col->getUnis(iTow);
    CalibData::TkrTotUni* dest;
    if (!(*unis)[uniIx]) {
      // have to make one
      (*unis)[uniIx]= makeUni(col);
    }
  
    dest = dynamic_cast<CalibData::TkrTotUni*>((*unis)[uniIx]);
    if (!dest) {
      // tilt!
      log << MSG::ERROR << "Unable to create CalibData::TkrTotUni object" 
          << endreq;
      return StatusCode::FAILURE;
    }  
    // For now, check if nstrips match.  If not, delete old 
    // (if exists) and allocate new.
    // But this isn't really quite right.  Can find out how
    // many strip objects are in ROOT object, but really what
    // we need to know is max strip id occurring.
    // Could add a "maxId"  field to calibRootData::TotUnilayer
    dest->resize(nStrips);

    if (nStrips > 0) {   // copy strip info
      for (unsigned iStrip = 0;  iStrip < (unsigned) nStrips; iStrip++) {
        const calibRootData::TotStrip* rStrip = rootUni->getStrip(iStrip);
        TkrTotStrip tdsStrip(rStrip->getStripId(), rStrip->getSlope(),
                             rStrip->getIntercept(), rStrip->getQuad(),
                             rStrip->getChi2(), rStrip->getDf());
        // bool ok = 
        dest->putStrip(tdsStrip);
      }
    }    
    delete rootUni;
    rootUni = 0;
  }    // end unilayer
  return StatusCode::SUCCESS;
}
