//$Header$
#include <stdio.h>
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/Service.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "CalibData/Moot/MootData.h"
#include "MootSvc/IMootSvc.h"
#include "CalibSvc/ICalibPathSvc.h"
#include <string>

/**
   @file UseMoot.cxx
   Simple algorithm to test functioning Moot data mediated by MootSvc
*/

  /** 
   @class UseMoot

   Algorithm exemplifying retrieval and use of Calorimeter pedestal calibration
*/
class UseMoot : public Algorithm {


public:
  UseMoot(const std::string& name, ISvcLocator* pSvcLocator); 

  StatusCode initialize();

  StatusCode execute();

  StatusCode finalize();

private:
  /// Helper function called by execute
  void processNewParm();
  void processNewFilters();

  IMootSvc*                m_pMootSvc;
  std::string              m_parmPath;
  unsigned                 m_hw;     // cached hardware key
  CalibData::MootParmCol*  m_pCol;
  bool                     m_latcParm;
  bool                     m_filterCfg;
  MsgStream*               m_log;
};


/// Instantiation of a static factory to create instances of this algorithm
//static const AlgFactory<UseMoot> Factory;
//const IAlgFactory& UseMootFactory = Factory;
DECLARE_ALGORITHM_FACTORY(UseMoot);

UseMoot::UseMoot(const std::string&  name, 
                 ISvcLocator*        pSvcLocator )
  : Algorithm(name, pSvcLocator), m_pMootSvc(0), 
    m_hw(0), m_pCol(0)
{
  // Declare properties here.
  Algorithm::declareProperty("FindLatcparm", m_latcParm = true);
  Algorithm::declareProperty("FindFilterCfg", m_filterCfg = true);
}


StatusCode UseMoot::initialize() {
  StatusCode sc;
  m_log = new MsgStream(msgSvc(), name());
  (*m_log) << MSG::INFO << "Initialize()" << endreq;

  // So far don't have any properties, but in case we do some day..
  setProperties();

  sc = service("MootSvc", m_pMootSvc, true);

  if ( !sc.isSuccess() ) {
    (*m_log) << MSG::ERROR 
             << "Could not get IMootSvc interface of MootSvc" 
             << endreq;
    return sc;
  }


  // Get properties from the JobOptionsSvc
  sc = setProperties();
  return StatusCode::SUCCESS;

}


StatusCode UseMoot::execute( ) {

  static bool firstCall = true;

  if (firstCall) {
    for (MOOT::InfoItem item = MOOT::INFOITEM_MOOTCONFIGKEY;
        item <= MOOT::INFOITEM_HWKEY; item = static_cast<MOOT::InfoItem> (item + 1) ) {
      MOOT::InfoSrc src = m_pMootSvc->getInfoItemSrc(item);
      (*m_log) << "Source for item " << item; 
      switch (src) {
      case MOOT::INFOSRC_TDS:
        (*m_log) << " is TDS";
        break;
      case MOOT::INFOSRC_JO:
        (*m_log) << " is job options";
        break;
      default:
        (*m_log) << " is unknown";
        break;
      }
      (*m_log) << endreq;
    }
    firstCall = false;
  }


  if (m_latcParm) processNewParm();
  if (m_filterCfg) processNewFilters();
  return StatusCode::SUCCESS;
}

void UseMoot::processNewFilters() {
  using CalibData::MootFilterCfg;

  static unsigned cfg = 0;

  (*m_log) << MSG::INFO << "Hi from UseMoot::processNewFilters" << std::endl;

  unsigned configKey = m_pMootSvc->getMootConfigKey();
  (*m_log) << MSG::INFO << "Using MOOT config #" << configKey << endreq;
  std::vector<MootFilterCfg> filters;
  if (cfg == configKey) return;
  else cfg = configKey;

  unsigned filterCnt = m_pMootSvc->getActiveFilters(filters);

  (*m_log) << MSG::INFO << "There are " << filterCnt << " active filters"
           << endreq;

  (*m_log) << MSG::INFO << "Filter names are: "  << endreq;
  for (unsigned ix = 0; ix < filterCnt; ix++) {
    (*m_log) << MSG::INFO << filters[ix].getName() << endreq;
  }

  filters.clear();

  unsigned acqMode = 1;

  filterCnt = m_pMootSvc->getActiveFilters(filters, acqMode);
  (*m_log) << MSG::INFO << "There are " << filterCnt 
           << " active filters for mode " << acqMode << endreq;

  for (acqMode = 0; acqMode < 8; acqMode ++) {
    for (unsigned hId = 0; hId < 8; hId++) {
      std::string hName;
      MootFilterCfg* pF = m_pMootSvc->getActiveFilter(acqMode, hId, hName);
      if (!pF) {
        (*m_log) << MSG::INFO << "No filter found for acqMode =" << acqMode
                 << " and handler id#" << hId << endreq;
      }
      else {
        (*m_log) << MSG::INFO << "For acqMode " << acqMode 
                 << " and handler id#" << hId << " handler name is: "
                 << hName << endreq
                 << "filter name is " << pF->getName() << endreq
                 << "schema id id is " << pF->getSchemaId() 
                 << " and instance id is " << pF->getInstanceId() << endreq;
      }
    }
  }
}
void UseMoot::processNewParm() {

  (*m_log) << MSG::INFO << "Hi from UseMoot::processNewParm" << std::endl;

  unsigned hw;

  std::string goodClass("latc_CFE_CAL_Mode"), badClass("latc_TEM");
  std::string goodPath, badPath;

  // Calling IMootSvc::getMootParmPath exercises essentially all
  // MootSvc client code
  goodPath = m_pMootSvc->getMootParmPath(goodClass, hw);
  if (hw == m_hw) {
    (*m_log) << "No change in hw key so nothing to do " << std::endl;
    return;
  }
  m_hw = hw;
  if (goodPath.size()) {
    (*m_log) << MSG::INFO << "Full path for parm class " << goodClass 
        << " is: " << goodPath << std::endl;
  }
  else {
    (*m_log) << MSG::INFO << "Full path for parm class " << goodClass
           << " not found " << std::endl;
  }
  badPath = m_pMootSvc->getMootParmPath(badClass, hw);
  if (badPath.size()) {
    (*m_log) << MSG::INFO << "Full path for parm class " << badClass 
        << " is: " << badPath << std::endl;
  }
  else (*m_log) << MSG::INFO << "Full path for parm class " << badClass
           << " not found " << std::endl;

  const CalibData::MootParm* gemParm = m_pMootSvc->getGemParm(hw);
  if (gemParm) {
    (*m_log) << "Gem parm src is " << gemParm->getSrc()
             << ", is from precinct " << gemParm->getPrecinct()
             << " and key is " << gemParm->getKey() << endreq;
  }
  else (*m_log) << "No gem parm found" << endreq;

  const CalibData::MootParm* roiParm = m_pMootSvc->getRoiParm(hw);
  if (roiParm) {
    (*m_log) << "Roi parm src is " << roiParm->getSrc()
             << ", is from precinct " << roiParm->getPrecinct()
             << " and key is " << roiParm->getKey() << endreq;
  }
  else (*m_log) << "No roi parm found" << endreq;
  
}

StatusCode UseMoot::finalize( ) {

  (*m_log) << MSG::INFO 
      << "          Finalize UseMoot "
      << endreq;
  
  return StatusCode::SUCCESS;
}

