// Include files
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/IDetDataSvc.h"

#include "CalibData/Cal/Gain.h"
#include "CalibData/Cal/Ped.h"
#include "CalibData/CalibTime.h"
#include "facilities/Timestamp.h"
 
#include "TestEnergyTool.h"

static ToolFactory<TestEnergyTool> s_factory;
const IToolFactory& TestEnergyToolFactory = s_factory;

TestEnergyTool::TestEnergyTool( const std::string& type, 
                                const std::string& name, 
                                const IInterface* parent)
  : AlgTool(type,name,parent) {
  declareInterface<ICalEnergyTool>(this);

  declareProperty( "startTime",  
                   m_startTimeAsc = "2003-1-10_00:20");
  declareProperty( "calibFlavor", m_calibFlavor = "none");
  
}

StatusCode TestEnergyTool::initialize() {
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "initialize" << endreq;

  StatusCode sc = StatusCode::SUCCESS;
        
  // extracting int constants
  double value;  // intermediate variable for reading constants from
  // xml file as doubles and converting them to interger 
  typedef std::map<int*,std::string> PARAMAP;
  PARAMAP param; // map containing pointers to integer constants to be read
  // with their symbolic names from xml file used as a key 

  //     filling the map with information on constants to be read 
    
  param[&m_pedestal]=std::string("cal.pedestal");
  param[&m_maxAdc]=std::string("cal.maxAdcValue");
  param[&m_thresh]=std::string("cal.zeroSuppressEnergy");
    
  // now try to find the GlastDevSvc service
    
  //    IGlastDetSvc* detSvc;
  sc = service("GlastDetSvc", detSvc);
    
  // loop over all constants information contained in the map
  for(PARAMAP::iterator it=param.begin(); it!=param.end();it++){

    //  attempt to get teh constant value via the method of GlastDetSvc
    if(!detSvc->getNumericConstByName((*it).second, &value)) {

      // if not successful - give the error message and return
      log << MSG::ERROR << " constant " <<(*it).second
          <<" not defined" << endreq;
      return StatusCode::FAILURE;

      //  if successful - fill the constant using the pointer from the map
    } else *((*it).first)= int(value);
  }
        
  // extracting double constants
    
  typedef std::map<double*,std::string> DPARAMAP;
  DPARAMAP dparam; // map containing pointers to double constants to be read
  // with their symbolic names from xml file used as a key 

  dparam[m_maxEnergy]=std::string("cal.maxResponse0");
  dparam[m_maxEnergy+1]=std::string("cal.maxResponse1");
  dparam[m_maxEnergy+2]=std::string("cal.maxResponse2");
  dparam[m_maxEnergy+3]=std::string("cal.maxResponse3");
    
  for(DPARAMAP::iterator dit=dparam.begin(); dit!=dparam.end();dit++){
    if(!detSvc->getNumericConstByName((*dit).second,(*dit).first)) {
      log << MSG::ERROR << " constant " <<(*dit).second << " not defined" << endreq;
      return StatusCode::FAILURE;
    } 
  }

	
  // get pointer to CalibDataSvc
  sc = service("CalibDataSvc", m_pCalibDataSvc, true);

  if ( !sc.isSuccess() ) {
    log << MSG::ERROR 
        << "Could not get IDataProviderSvc interface of CalibDataSvc" 
        << endreq;
    return sc;
  }
 
  // Query the IDetDataSvc interface of the calib data service
  sc = m_pCalibDataSvc->queryInterface(IID_IDetDataSvc, 
                                       (void**) &m_detDataSvc);
  if ( !sc.isSuccess() ) {
    log << MSG::ERROR 
        << "Could not query IDetDataSvc interface of CalibDataSvc" 
        << endreq;
    return sc;
  } else {
    log << MSG::DEBUG 
        << "Retrieved IDetDataSvc interface of CalibDataSvc" 
        << endreq;
  }
        
  // Get properties from the JobOptionsSvc
  sc = setProperties();
  if ( !sc.isSuccess() ) {
    log << MSG::ERROR << "Could not set jobOptions properties" << endreq;
    return sc;
  }

  std::string strTimeAsc = m_startTimeAsc.value();
  unsigned int underpos = strTimeAsc.find("_");
  if (underpos < strTimeAsc.size())
    strTimeAsc.replace(underpos, 1, " ");
  
  m_startTime = facilities::Timestamp(strTimeAsc).getClibTime();

  log << MSG::DEBUG << "Properties were read from jobOptions" << endreq;
  log << MSG::INFO << "Time of first event: (ascii) "
      << strTimeAsc       << endreq; 
  log << MSG::INFO << "Time of first event: (seconds since 1970) "
      << m_startTime       << endreq; 

  //return sc;

 
  return StatusCode::SUCCESS;
}

/*!
  Method:
  -# if(calibFlavor != none) set event time for DetDataSvc, obtain pointers pedestal and elecgain data in TDS, use uniqe pedestal and gain information for each xtal/range
  -# else use default pedesatl values for all xtal/range</ol>
  -# convert adc values into energy (e = gain*(adc-ped))
  -# return average of POS and NEG energy estimations

  \note The replication of the original CalXtalRec algorithm is incomplete as ICalEnergyTool interface does not
  allow for returning two separate energy values.  In this case the average energy (POS+NEG/2) is returned. Those who 
  need energy from both faces may want to call the overloaded version of the function which deals with each
  face individually.

*/
StatusCode TestEnergyTool::calculate(const idents::CalXtalId &xtalId, 
                                     idents::CalXtalId::AdcRange rangeP,
                                     idents::CalXtalId::AdcRange rangeN,
                                     int adcP, 
                                     int adcN,
                                     float position,
                                     float &energy                    // output
                                     ) {
  MsgStream log(msgSvc(), name());
  pPeds = 0;
  pGains = 0;

  if(m_calibFlavor != "none") {

    // Set the event time
    facilities::Timestamp time = facilities::Timestamp(m_startTime);
    log << MSG::DEBUG << "Event time: "
        << time.getString()
        << endreq; 
    CalibData::CalibTime ctime(time);
    log << MSG::DEBUG << "Event time (hours) " << ctime.hours() << endreq;
    m_detDataSvc->setEventTime(ctime);
  
    std::string fullPedPath = "/Calib/CAL_Ped/"+m_calibFlavor;
    std::string fullGainPath = "/Calib/CAL_ElecGain/"+m_calibFlavor;

    DataObject *pObject;
        
    //getting pointers to the calibration data of each type
    if(m_pCalibDataSvc->retrieveObject(fullPedPath, pObject) == StatusCode::SUCCESS) {

      pPeds = dynamic_cast<CalibData::CalCalibPed *> (pObject);
      if (!pPeds) {
        log << MSG::ERROR << "Dynamic cast to CalCalibPed failed" << endreq;
        return StatusCode::FAILURE;
      }
    } else{
      log << MSG::INFO << "Enable to retrieve pedestals from calib database" << endreq;
    }
    if(m_pCalibDataSvc->retrieveObject(fullGainPath, pObject) == StatusCode::SUCCESS) {

      pGains = dynamic_cast<CalibData::CalCalibGain *> (pObject);
      if (!pGains) {
        log << MSG::ERROR << "Dynamic cast to CalCalibGain failed" << endreq;
        return StatusCode::FAILURE;
      }
    }else{            
      log << MSG::INFO << "Enable to retrieve gains from calib database" << endreq;
    }

  }

  //extraction of pedestals
  
  // default value of pedestal - if calibration database isn't accessible
  float pedP = m_pedestal;
  float pedN = m_pedestal;
  
  if(pPeds){
    CalibData::RangeBase* pRangeP = 
      pPeds->getRange(xtalId, rangeP,idents::CalXtalId::POS);
    CalibData::RangeBase* prangeN = 
      pPeds->getRange(xtalId, rangeN,idents::CalXtalId::NEG);
    
    CalibData::Ped* pPedP = 
      dynamic_cast<CalibData::Ped * >(pRangeP);
    CalibData::Ped* ppedN = 
      dynamic_cast<CalibData::Ped * >(prangeN);
    pedP = pPedP->getAvr();
    pedN = ppedN->getAvr();
  }
  
  // extraction of gains
  // default values of gains -if calibration database isn't accessible
  float gainP = m_maxEnergy[rangeP]/(m_maxAdc-m_pedestal);
  float gainM = m_maxEnergy[rangeN]/(m_maxAdc-m_pedestal);
  
  if(pGains){
    CalibData::RangeBase* pRangeP = pGains->getRange(xtalId, rangeP,idents::CalXtalId::POS);
    CalibData::RangeBase* prangeN = pGains->getRange(xtalId, rangeN,idents::CalXtalId::NEG);
    
    CalibData::Gain* pGainP = dynamic_cast<CalibData::Gain * >(pRangeP);
    CalibData::Gain* pGainM = dynamic_cast<CalibData::Gain * >(prangeN);
      
    gainP = pGainP->getGain();
    gainM = pGainM->getGain();
  }
  
  // convert adc values into energy
  double eneP = gainP*(adcP-pedP);
  double eneM = gainM*(adcN-pedN);
  
  log << MSG::INFO << " eneP " << eneP << " eneM " << eneM << endreq;
  energy = (eneP+eneM)/2;

  return StatusCode::SUCCESS;
}
  
/*!
  Method:
  -# if(calibFlavor != none) set event time for DetDataSvc, obtain pointers pedestal and elecgain data in TDS, use uniqe pedestal and gain information for each xtal/range
  -# else use default pedestal values for all xtal/range
  -# convert adc values into energy (e = gain*(adc-ped))
*/
StatusCode TestEnergyTool::calculate(const idents::CalXtalId &xtalId,
                                     idents::CalXtalId::AdcRange range,
                                     idents::CalXtalId::XtalFace face,
                                     int adc, 
                                     float position,
                                     float &energy                    // output
                                     ) {

  MsgStream log(msgSvc(), name());
  pPeds = 0;
  pGains = 0;

  if(m_calibFlavor != "none") {
    // Set the event time
    facilities::Timestamp time = facilities::Timestamp(m_startTime);
    log << MSG::DEBUG << "Event time: "
        << time.getString()
        << endreq; 
    CalibData::CalibTime ctime(time);
    log << MSG::DEBUG << "Event time (hours) " << ctime.hours() << endreq;
    m_detDataSvc->setEventTime(ctime);
  
    std::string fullPedPath = "/Calib/CAL_Ped/"+m_calibFlavor;
    std::string fullGainPath = "/Calib/CAL_ElecGain/"+m_calibFlavor;

    DataObject *pObject;
        
    //getting pointers to the calibration data of each type
    if(m_pCalibDataSvc->retrieveObject(fullPedPath, pObject) == StatusCode::SUCCESS) {

      pPeds = dynamic_cast<CalibData::CalCalibPed *> (pObject);
      if (!pPeds) {
        log << MSG::ERROR << "Dynamic cast to CalCalibPed failed" << endreq;
        return StatusCode::FAILURE;
      }
    } else{
      log << MSG::INFO << "Enable to retrieve pedestals from calib database" << endreq;
    }
    if(m_pCalibDataSvc->retrieveObject(fullGainPath, pObject) == StatusCode::SUCCESS) {

      pGains = dynamic_cast<CalibData::CalCalibGain *> (pObject);
      if (!pGains) {
        log << MSG::ERROR << "Dynamic cast to CalCalibGain failed" << endreq;
        return StatusCode::FAILURE;
      }
    }else{            
      log << MSG::INFO << "Enable to retrieve gains from calib database" << endreq;
    }

  }

  //extraction of pedestals
  
  // default value of pedestal - if calibration database isn't accessible
  float ped = m_pedestal;
  
  if(pPeds){
    CalibData::RangeBase* pRange = 
      pPeds->getRange(xtalId, range,face);
    
    CalibData::Ped* pPed = 
      dynamic_cast<CalibData::Ped * >(pRange);
    ped = pPed->getAvr();
  }
  
  
  // extraction of gains
  // default values of gains -if calibration database isn't accessible
  float gain = m_maxEnergy[range]/(m_maxAdc-m_pedestal);
    
  if(pGains){
    CalibData::RangeBase* pRange = pGains->getRange(xtalId, range,face);
        
    CalibData::Gain* pGain = dynamic_cast<CalibData::Gain * >(pRange);
          
    gain = pGain->getGain();
  }
  
  
  // convert adc values into energy
  double ene = gain*(adc-ped);
    
  energy = ene;

  return StatusCode::SUCCESS;
}
