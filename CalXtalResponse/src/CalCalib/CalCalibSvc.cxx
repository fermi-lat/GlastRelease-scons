// $Header $
/** @file
    @author Z.Fewtrell
 */
// LOCAL
#include "CalCalibSvc.h"

// GLAST

// EXTLIB
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SvcFactory.h"
#include "GaudiKernel/Incident.h"
#include "GaudiKernel/IIncidentSvc.h"
#include "Event/TopLevel/Event.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "facilities/Timestamp.h"

#include <cmath>

// STD
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <functional>

using namespace CalUtil;

//static SvcFactory< CalCalibSvc > a_factory;
//const ISvcFactory& CalCalibSvcFactory = a_factory; 
DECLARE_SERVICE_FACTORY(CalCalibSvc);

CalCalibSvc::CalCalibSvc(const string& name, ISvcLocator* Svc) 
  : Service(name,Svc),
    m_nEvent(0),
    m_pedT0(CalUtil::RngIdx::N_VALS,-5000),
    m_pedTempCoef(CalUtil::RngIdx::N_VALS,0),
    m_pedMgr(m_ccsShared),
    m_inlMgr(m_ccsShared),
    m_asymMgr(m_ccsShared),
    m_mpdMgr(m_ccsShared),
    m_tholdCIMgr(m_ccsShared),
    m_eventSvc(0),
    m_gotPedTime(false)
{
  // declare the properties
  declareProperty("CalibDataSvc",      m_ccsShared.m_calibDataSvcName = 
                  "CalibDataSvc");
  declareProperty("idealCalibXMLPath", m_ccsShared.m_idealCalibXMLPath = 
                  "$(CALUTILXMLPATH)/idealCalib_flight.xml");
  declareProperty("DefaultFlavor", m_defaultFlavor    
                  = "ideal");
  declareProperty("FlavorIntNonlin", m_flavorIntNonlin  = "");
  declareProperty("FlavorAsym",      m_flavorAsym       = "");
  declareProperty("FlavorPed",       m_flavorPed        = "");
  declareProperty("FlavorMeVPerDac", m_flavorMPD        = "");
  declareProperty("FlavorTholdCI",   m_flavorTholdCI    = "");
  declareProperty("TemperatureFile",   m_temperatureFile    = "");
  declareProperty("PedTempCor",   m_pedTempCorFile    = "");
}

StatusCode  CalCalibSvc::queryInterface (const InterfaceID& riid, void **ppvIF) {
  if (IID_ICalCalibSvc == riid) {
    *ppvIF = (ICalCalibSvc*)(this);
    return StatusCode::SUCCESS;
  } else return Service::queryInterface (riid, ppvIF);
}

/// intialize / retrieve all needed Gaudi based objects
StatusCode CalCalibSvc::initialize () 
{
  // Call super-class
  Service::initialize ();

  MsgStream msglog(msgSvc(), name()); 
  msglog << MSG::INFO << "initialize" << endreq;

  StatusCode sc;

  //-- jobOptions --//
  if ((sc = setProperties()).isFailure()) {
    msglog << MSG::ERROR << "Failed to set properties" << endreq;
    return sc;
  }

  // Post-init processing of jobOptions.
  // set default flavors unless otherwise specified
  if (!m_flavorIntNonlin.value().length())  m_flavorIntNonlin = m_defaultFlavor;
  if (!m_flavorAsym.value().length())       m_flavorAsym      = m_defaultFlavor;
  if (!m_flavorPed.value().length())        m_flavorPed       = m_defaultFlavor;
  if (!m_flavorMPD.value().length())        m_flavorMPD       = m_defaultFlavor;
  if (!m_flavorTholdCI.value().length())    m_flavorTholdCI   = m_defaultFlavor;

  msglog << MSG::DEBUG << "Initializing..."     << endreq;
  msglog << MSG::DEBUG << "  CalibDavaSvc   : " << m_ccsShared.m_calibDataSvcName  << endreq;
  msglog << MSG::DEBUG << "  DefaultFlavor  : " << m_defaultFlavor     << endreq;
  msglog << MSG::DEBUG << "  FlavorIntNonlin: " << m_flavorIntNonlin   << endreq;
  msglog << MSG::DEBUG << "  FlavorAsym     : " << m_flavorAsym        << endreq;
  msglog << MSG::DEBUG << "  FlavorPed      : " << m_flavorPed         << endreq;
  msglog << MSG::DEBUG << "  FlavorMeVPerDac: " << m_flavorMPD         << endreq;  
  msglog << MSG::DEBUG << "  FlavorTholdCI  : " << m_flavorTholdCI     << endreq;    

  msglog << MSG::INFO << "  TemperatureFile  : " << m_temperatureFile     << endreq;    
  msglog << MSG::INFO << "  PedTempCor  : " << m_pedTempCorFile     << endreq;    


  /// init all the CalCalibShared stuff
  sc = m_ccsShared.initialize(*this);
  if (sc.isFailure()) return sc;

  // Initialize individual CalibItemMgr members.
  sc = m_mpdMgr.initialize(m_flavorMPD);
  if (sc.isFailure()) return sc;
  sc = m_pedMgr.initialize(m_flavorPed);
  if (sc.isFailure()) return sc;
  sc = m_asymMgr.initialize(m_flavorAsym);
  if (sc.isFailure()) return sc;
  sc = m_inlMgr.initialize(m_flavorIntNonlin);
  if (sc.isFailure()) return sc;
  sc = m_tholdCIMgr.initialize(m_flavorTholdCI);
  if (sc.isFailure()) return sc;


  if(m_temperatureFile.value() != ""){

    ifstream tempfile(m_temperatureFile.value().c_str());
    if (!tempfile.is_open())
            throw runtime_error(string("Unable to open " + m_temperatureFile.value()));
    string line;
    while (tempfile.good()) {
      unsigned short day;
      unsigned short month;
      unsigned short year;
      unsigned short hour;
      unsigned short minutes;
      unsigned short seconds;
      getline(tempfile, line);
      if (tempfile.fail()) break; // bad get
    
      // check for comments
      if (line[0] == ';')
        continue;

      istringstream istrm(line);
      
      istrm >> day >> month >> year >> hour >> minutes >> seconds;

      int time = ((day+3 + (month-8)*31)*24 + hour)*60+minutes;
      m_tempTime.push_back(time);
      for (int twr=3;twr>0;twr--){
	float temp0,temp1,temp2,temp3,tempav;
	istrm >> temp0 >> temp1 >> temp2 >> temp3;
	tempav=(temp0+temp1+temp2+temp3)/4.0;
	float dtempSPS[]={0,1.0,1.3,0.7};
	if(time/1440 > 30)tempav -= dtempSPS[twr];  //   temperature bias during SPS beam test period
	(m_cuTwrTemp[twr]).push_back(tempav);
      }
    }
    msglog << MSG::INFO << "  successfully read " << m_tempTime.size()
	   << " temperature records "     << endreq;   
  }
 
  if(m_pedTempCorFile.value() != ""){

    ifstream pedfile(m_pedTempCorFile.value().c_str());
    if (!pedfile.is_open())
            throw runtime_error(string("Unable to open " + m_pedTempCorFile.value()));
    string line;
    while (pedfile.good()) {
      unsigned short twr;
      unsigned short lyr;
      unsigned short col;
      unsigned short face;
      unsigned short rng;
      float ped,dpdt;
      
      getline(pedfile, line);
      if (pedfile.fail()) break; // bad get

      // check for comments
      if (line[0] == ';')
        continue;

      istringstream istrm(line);

      istrm >> twr
            >> lyr
            >> col
            >> face
            >> rng
            >> ped
            >> dpdt;

      //      cout <<" "<< twr <<" "  << lyr <<" " << col <<" "
      //	   << face <<" " << rng<<" " << ped<<" " << dpdt << endl;

      const RngIdx rngIdx(twr,
                    LyrNum(lyr),
                    col,
                    FaceNum((idents::CalXtalId::XtalFace)face),
                    RngNum(rng));

      m_pedT0[rngIdx]   = ped;
      m_pedTempCoef[rngIdx] = dpdt;

    }

    msglog << MSG::INFO << "  successfully read pedestal temperature corrections"      
	   << endreq;   

  }


  // Get ready to listen for BeginEvent
  IIncidentSvc* incSvc;
  sc = service("IncidentSvc", incSvc, true);
  if (sc.isSuccess() ) {
    incSvc->addListener(this, "BeginEvent", INCIDENT_PRIORITY);
  } else {
    msglog << MSG::ERROR << "can't find IncidentSvc" << endreq;
    return sc;
  }


// Need event data service for timestamp stuff
  sc = serviceLocator()->service("EventDataSvc", m_eventSvc, true);
  if (sc .isFailure() ) {
    msglog << MSG::ERROR << "Unable to find EventDataSvc " << endreq;
    return sc;
  }


  return StatusCode::SUCCESS;
}

/// Inform child objects that a new event has occured
void CalCalibSvc::handle ( const Incident& inc ) { 
  if ((inc.type() == "BeginEvent")) {
    m_mpdMgr.invalidate();
    m_pedMgr.invalidate();
    m_asymMgr.invalidate();
    m_inlMgr.invalidate();
    m_tholdCIMgr.invalidate();

    m_gotPedTime = false;
  }
  return; 
}

StatusCode CalCalibSvc::evalFaceSignal(const RngIdx rngIdx, const float adcPed, float &ene) {
  StatusCode sc;

  // adc -> cidac
  float cidac;
  sc = evalCIDAC(rngIdx, adcPed, cidac);
  if (sc.isFailure()) return sc;

  float mpdDiode;
  sc = getMPDDiode(rngIdx.getDiodeIdx(), mpdDiode);
  if (sc.isFailure()) return StatusCode::FAILURE;

  ene = cidac*mpdDiode;
  
  return StatusCode::SUCCESS;
}


  StatusCode CalCalibSvc::getPed(CalUtil::RngIdx rngIdx,float& ped)
{

  static const facilities::Timestamp bt06Start("2006-7-28 00:00");
  static const facilities::Timestamp missionStart("2001-1-1 00:00");
  static const unsigned bt06Sec = (unsigned) bt06Start.getClibTime();
  static const unsigned missionSec = (unsigned) missionStart.getClibTime();

  MsgStream msglog(msgSvc(), name()); 

  if (m_gotPedTime == false) {
    SmartDataPtr<Event::EventHeader> eventHeader(m_eventSvc, "/Event");
  
    if (!eventHeader) {
      msglog << MSG::ERROR << "Unable to retrieve event timestamp for digis" 
             << endreq;
      return 0;
    }
    
    m_nEvent++;
 
  if(m_pedTempCorFile.value() != "" && m_temperatureFile.value() != ""){

    double fromMissionStart = (eventHeader->time()).time();
    int fromBtStartSec = fromMissionStart-(bt06Sec-missionSec);
    vector<int>::iterator it_btstart = m_tempTime.begin();
    vector<int>::iterator it_btend = m_tempTime.end();
    vector<int>::iterator it_curtime = lower_bound(it_btstart,it_btend, fromBtStartSec/60);
    unsigned cur_temp_idx = it_curtime - it_btstart;
    for(int twr=1;twr<4;twr++)
      m_cur_temp_twr[twr] = (m_cuTwrTemp[twr])[cur_temp_idx];



    if ((m_nEvent == ((m_nEvent/1000) * 1000) ) ) {

      msglog << MSG::DEBUG << "event time in seconds from beam test start: " 
	     << fromBtStartSec 
                      <<" temp_twr1=" << m_cur_temp_twr[1]
	              <<" temp_twr2=" << m_cur_temp_twr[2] 
	              <<" temp_twr3=" << m_cur_temp_twr[3]
	     << endreq;
    }
  }
    m_gotPedTime = true;
  }

  const CalibData::Ped* pedptr = m_pedMgr.getPed(rngIdx);
  if(pedptr == 0)return StatusCode::FAILURE;
	       
  ped = pedptr->getAvr();

  if(m_pedTempCorFile.value() != ""){

    float ped0 = m_pedT0[rngIdx];
    if(ped0 > 0){
      float dpdt = m_pedTempCoef[rngIdx];
      int twr = rngIdx.getTwr().val();
      if(twr>=0 && twr<=3) ped = ped0+dpdt*(m_cur_temp_twr[twr]-23.0);
    }
  }
  return StatusCode::SUCCESS;
}

StatusCode CalCalibSvc::getPedSig(CalUtil::RngIdx rngIdx, float &sig){

  const CalibData::Ped* pedptr = m_pedMgr.getPed(rngIdx);
  if(pedptr == 0)return StatusCode::FAILURE;
	       
  sig = pedptr->getSig();
  return StatusCode::SUCCESS;
}



StatusCode CalCalibSvc::getMPDDiode(const CalUtil::DiodeIdx diodeIdx, float &mpdDiode) {
  const XtalIdx xtalIdx(diodeIdx.getXtalIdx());

  // MeVPerDAC
  // need to create tmp rngIdx w/ only twr/lyr/col info
  CalibData::CalMevPerDac const*const calMPD = getMPD(xtalIdx);
  if (!calMPD) return StatusCode::FAILURE;

  float mpd;
  float asymCtr;

  // get overall asymmetry & mpd for this xtal
  StatusCode sc;
  if (diodeIdx.getDiode() == LRG_DIODE) {
    mpd = calMPD->getBig()->getVal();
    sc = getAsymCtr(xtalIdx, ASYM_LL, asymCtr);
    if (sc.isFailure()) return sc;
  }
  else { // diode == SM_DIODE
    mpd = calMPD->getSmall()->getVal();
    sc = getAsymCtr(xtalIdx, ASYM_SS, asymCtr);
    if (sc.isFailure()) return sc;
  }

  // correct for overall asymmetry of diodes (use asym at center
  // of xtal)
  if (diodeIdx.getFace() == POS_FACE)
    mpd *= exp(-1*asymCtr/2);
  else // face == NEG_FACE
    mpd *= exp(asymCtr/2);

  mpdDiode = mpd;


  return StatusCode::SUCCESS;
}



