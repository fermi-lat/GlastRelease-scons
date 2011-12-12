// $Header$
/**
            @file  RootBaseCnv.cxx

   Implementation file for Root calibration converter base class
*/

#include "RootBaseCnv.h"

#include "GaudiKernel/CnvFactory.h"
#include "GaudiKernel/IOpaqueAddress.h"
#include "GaudiKernel/DataObject.h"
#include "GaudiKernel/IAddressCreator.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/IConversionSvc.h"
#include "GaudiKernel/MsgStream.h"

#include "facilities/Util.h"   // for translating env variables
#include "CalibSvc/ICalibRootSvc.h" 
#include "CalibSvc/ICalibMetaCnvSvc.h"
#include "CalibSvc/IInstrumentName.h"
#include "CalibData/CalibBase.h"

// Guessing at needed Root includes
#include "TROOT.h"    // need this for cd??
#include "TFile.h"
#include "TTree.h"
#include "TObject.h"



// May need the following some day for dac collections
// #include "commonRootData/idents/CalXtalId.h"
// #include CalibData/DacCol.h"


RootBaseCnv::~RootBaseCnv() {
  // release TFile, TTree if they need releasing.  With normal 
  // termination they should already have been released.

  //  doClean();

}

// static CnvFactory<RootBaseCnv> s_factory;
// const ICnvFactory& RootBaseCnvFactory = s_factory;
RootBaseCnv::RootBaseCnv( ISvcLocator* svc, const CLID& clid) :
  Converter (CALIBROOT_StorageType, clid, svc),
  m_rootSvc (0), m_metaSvc(0), m_instrSvc(0), m_vstart(0), m_vend(0),
  m_outFile(0), m_ttree(0), m_inFile(0), m_saveDir(0) {}

StatusCode RootBaseCnv::initialize() {
  StatusCode status = Converter::initialize();

  IDataProviderSvc* dp;

  // I guess the service names are assigned in jobOptions?

  serviceLocator()->getService ("CalibDataSvc",
                                IDataProviderSvc::interfaceID(),
                                (IInterface*&)dp);
  setDataProvider(dp);
  
  // Locate the Root Conversion Service
  serviceLocator()->getService ("CalibRootCnvSvc",
                                IID_ICalibRootSvc,
                                (IInterface*&) m_rootSvc);

  // Locate meta conversion service
  // Will anything need to be changed here to accommodate possibility
  // of two concrete implementations of ICalibMetaCnvSvc?  Would
  // have different storage types.  Could specify type desired
  // as job option.  Ditto for name of class?
  serviceLocator()->getService("CalibMySQLCnvSvc", 
                               IID_ICalibMetaCnvSvc,
                                (IInterface*&)m_metaSvc);

  serviceLocator()->getService ("CalibDataSvc",
                                IID_IInstrumentName,
                                (IInterface*&)m_instrSvc);

  return status;
}

StatusCode RootBaseCnv::finalize() {
  return Converter::finalize();
}


            /******      ROOT services     *****/

StatusCode RootBaseCnv::createRoot(const std::string& /* fname */, 
                                   CalibData::CalibBase* /* pTDSObj */) {
  MsgStream log(msgSvc(), "RootBaseCnv");
  log << MSG::ERROR 
      << "createRoot method not implemented for this calibration type" 
      << endreq;
  return StatusCode::FAILURE;
}

/*
StatusCode RootBaseCnv::openRead(const std::string& fname, 
                                 const std::string& branch,
                                 TObject*& pCalib) {
*/
StatusCode RootBaseCnv::openRead(const std::string& fname) { 

  MsgStream log(msgSvc(), "RootBaseCnv");

  // Check fname isn't empty
  if (fname == std::string("")) return StatusCode::FAILURE;

  if (doClean() ) {
    log << MSG::WARNING << "Previous operation didn't clean up! " << endreq;
  }
  m_saveDir = gDirectory;

  std::string ourName(fname);
  facilities::Util::expandEnvVar(&ourName);

  //  TFile  f(ourName.c_str());
  m_inFile = new TFile(ourName.c_str());
  
  if (!m_inFile->IsOpen() ) {
    log << MSG::ERROR << "ROOT file " << ourName 
        << "could not be opened for reading " << endreq;
    delete m_inFile;
    m_inFile = 0;
    // if (m_crash)  {   don't have an m_crash
    log << MSG::FATAL << std::endl << "Exiting... " << std::endl << endreq;
    exit(1);
    // }

    return StatusCode::FAILURE;
  }
  else {
    log << MSG::INFO
        << "Successfully opened ROOT file " << fname << " aka " << ourName
        << " for reading " << endreq;
  }


  m_inFile->cd();             //    Maybe will need this

  //  ##### This may only be appropriate for CAL files ####
  //  ##### Probably should do it in RootBaseCalCnv
  /*
  TTree* pTree = (TTree*)m_inFile->Get("Calib");

  pTree->SetBranchAddress(branch.c_str(), &pCalib);
  pTree->GetEvent();
 
  //  Is this really what I need to do?  Or just a shortcut
  // while I concentrate on Write?
  // Maybe it is -- at this point information has been read from file
  // into TObject.
  m_inFile->Close();
  m_saveDir->cd();
  */
  return StatusCode::SUCCESS;
}

StatusCode RootBaseCnv::closeRead() {
  m_inFile->Close();

  delete m_inFile;
  m_inFile = 0;

  if (m_saveDir) {
    m_saveDir->cd();
    m_saveDir = 0;
  }
  return StatusCode::SUCCESS;
}
  
StatusCode RootBaseCnv::openWrite(const std::string& fname, 
                                  const std::string& className,
                                  TObject*& pCalib) {

  MsgStream log(msgSvc(), "RootBaseCnv");

  // Check fname isn't empty
  if (fname == std::string("")) return StatusCode::FAILURE;

  std::string ourName(fname);
  facilities::Util::expandEnvVar(&ourName);

  if (doClean() ) {
    log << MSG::WARNING << "Previous operation didn't clean up! " << endreq;
  }

  m_saveDir = gDirectory;

  // Should also check that there is no in-progress read.  OR else
  // make it possible to have them both going
  m_outFile = new TFile(ourName.c_str(), "RECREATE");
  if (!m_outFile->IsOpen()) {
    log << MSG::ERROR << "ROOT file " << fname << " aka " << ourName 
        << " could not be opened for writing" << endreq;
    delete m_outFile;
    m_outFile = 0;
    return StatusCode::FAILURE;
  }
  else {
    log << MSG::INFO
        << "Successfully opened ROOT file " << fname << " aka " << ourName
        << " for writing " << endreq;
  }
  m_outFile->cd();

  // Our tree is always called "Calib".  Branch is named after calib. type
  m_ttree = new TTree("Calib", "GLAST calibration data");

  m_ttree->Branch(className.c_str(), className.c_str(), &pCalib);
  return StatusCode::SUCCESS;
}

StatusCode RootBaseCnv::closeWrite() {

  MsgStream log(msgSvc(), "RootBaseCnv");

  StatusCode ret = StatusCode::SUCCESS;

  m_ttree->Fill();
  m_outFile->cd();
  int nBytes = m_outFile->Write(0, TObject::kOverwrite);
  if (!nBytes) {
    log << MSG::ERROR << "Unable to write file " << endreq;
    ret = StatusCode::FAILURE;
  }
  m_outFile->Close();
  delete m_outFile;
  m_outFile = 0;
  if (m_saveDir) m_saveDir->cd();
  m_saveDir = 0;
  return ret;
}


StatusCode RootBaseCnv::readRootObj(const std::string& treename, 
                                    const std::string& branch,
                                    TObject*& pObj, unsigned ix){
  TTree* pTree = (TTree*)m_inFile->Get(treename.c_str());

  return readRootObj(pTree, branch, pObj, ix);
 }

StatusCode RootBaseCnv::readRootObj(TTree* pTree,
                                    const std::string& branch,
                                    TObject*& pObj, unsigned ix){
  TBranch* pBranch=pTree->GetBranch(branch.c_str());
  pBranch->SetAddress(&pObj);
  int nBytes = pBranch->GetEntry(ix);
  return (nBytes > 0) ? StatusCode::SUCCESS : StatusCode::FAILURE;
 }

bool RootBaseCnv::doClean() {
  bool ret = false;

  if (m_outFile) {
    m_outFile->Close();
    delete m_outFile;
    m_outFile = 0;
    ret = true;
  }
  if (m_inFile) {
    m_inFile->Close();
    delete m_inFile;
    m_inFile = 0;
    ret = true;
  }
  m_ttree = 0;
  if (m_saveDir) {
    m_saveDir->cd();
    m_saveDir = 0;
  }
  return ret;
}


// Do our part to write out object -- which is nothing
StatusCode RootBaseCnv::fillRoot(CalibData::CalibBase* /*pTDSObj */, 
                                 TObject* /* pRootObj */) {

  // Get instrument name from InstrumentName service  Now handled by 
  // RootCalBaseCnv
  //  TString instr = TString((m_instrSvc->getInstrumentName()).c_str());
  //  pRootObj->setInstrument(instr);
  return StatusCode::SUCCESS;
}

// (To TDS) Conversion stuff
StatusCode RootBaseCnv::createObj(IOpaqueAddress* addr,
                                  DataObject*& refpObject) {
  //  StatusCode ret;

  // first do the things we always need:
  //   First string parameter of opaque address is file ident
  const std::string* par = addr->par();

  std::string par0 = par[0];

  return internalCreateObj(par0, refpObject, addr);

}

StatusCode RootBaseCnv::internalCreateObj(const std::string& fname,
                                          DataObject*& refpObject,
                                          IOpaqueAddress* address) {
  MsgStream log(msgSvc(), "RootBaseCnv");
  RootBaseCnv* converter = this;
  CLID classId = address->clID();
  IConverter* conv = this->conversionSvc()->converter(classId);

  if (0 == conv) {
    log << MSG::WARNING
        << "No proper converter found for classID " << classId
            << ", the default converter"
            << " will be used. " << endreq;
  } else {
    converter = dynamic_cast <RootBaseCnv*> (conv);
    if (0 == converter) {
      log << MSG::ERROR
          << "The converter found for classID " << classId
              << " was not a descendent of RootBaseCnv as it should be "
              << "( was of type " << typeid (*converter).name() << "). "
              << "The default converter will be used" << endreq;
      converter = this;
    }
  }

  unsigned int serNo = *(address->ipar());
  m_serNo = serNo;
  StatusCode sc = m_metaSvc->getValidInterval(serNo, 
                                              &m_vstart, 
                                              &m_vend );

  // creates an object for the node found
  if (sc.isSuccess()) sc = converter->i_createObj(fname, refpObject);
  if (sc.isFailure()) {
    return sc;
  }

  // ends up the object construction
  sc = converter->i_processObj(refpObject, address);
  if (sc.isSuccess()) {
    log << MSG::DEBUG << "Successfully created calib. object " << endreq;
  }
  return sc;
}

/* 
   Base class version of this routine shouldn't really be called
   since it doesn't correspond to a TDS object.
*/
StatusCode RootBaseCnv::i_createObj (const std::string& /* fname */,
                                     DataObject*& /* refpObject */) {
  return StatusCode::FAILURE;  // shouldn't ever get here
}
  
// Default is to do nothing.  Derived classes may override.
StatusCode RootBaseCnv::i_processObj(DataObject*, // pObject,
                                     IOpaqueAddress* ) /* address */  {
  return StatusCode::SUCCESS;
}

/// Another utility for derived classes to use
void RootBaseCnv::setBaseInfo(CalibData::CalibBase* pObj) {
  pObj->setValidity(*m_vstart, *m_vend);
  pObj->setSerNo(m_serNo);
}
