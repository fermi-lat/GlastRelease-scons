//    $Header$

/** @file     implement XtalDigiTool.h
    @author Zach Fewtrell

*/

// Include files

// LOCAL
#include "XtalDigiTool.h"
#include "../CalFailureMode/CalFailureModeSvc.h"
#include "CalXtalResponse/ICalCalibSvc.h"
#include "../CalCalib/IPrecalcCalibTool.h"

// GLAST
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "idents/VolumeIdentifier.h"
#include "Event/TopLevel/Event.h"
#include "Event/TopLevel/EventModel.h"
#include "Event/Digi/CalDigi.h"


// EXTLIB
#include "GaudiKernel/ToolFactory.h"
#include "TTree.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/IDataProviderSvc.h"

// STD
#include <map>

using namespace CalUtil;
using namespace Event;
using namespace CalibData;
using namespace idents;
using namespace std;

static ToolFactory<XtalDigiTool> s_factory;
const IToolFactory& XtalDigiToolFactory = s_factory;

static float round_int(float in) { return floor(in + 0.5);}

XtalDigiTool::XtalDigiTool( const string& type,
                            const string& name,
                            const IInterface* parent)
  : AlgTool(type,name,parent),
    m_calCalibSvc(0),
    m_maxAdc(-1),
    m_tuple(0),
    m_calFailureModeSvc(0),
    m_evtSvc(0),
    m_precalcCalib(0)
{
  declareInterface<IXtalDigiTool>(this);

  declareProperty("CalCalibSvc",        m_calCalibSvcName    = "CalCalibSvc");
  declareProperty("PrecalcCalibTool",   m_precalcCalibName   = "PrecalcCalibTool");
  declareProperty("tupleFilename",      m_tupleFilename      = "");
  declareProperty("tupleLACOnly",       m_tupleLACOnly       = true);

}

StatusCode XtalDigiTool::initialize() {
  MsgStream msglog(msgSvc(), name());   
  msglog << MSG::INFO << "initialize" << endreq;

  StatusCode sc;

  double val;
  typedef std::map<int*,string> PARAMAP;

  //-- jobOptions --//
  if ((sc = setProperties()).isFailure()) {
    msglog << MSG::ERROR << "Failed to set properties" << endreq;
    return sc;
  }

  // try to find the GlastDevSvc service
  IGlastDetSvc* detSvc;
  sc = service("GlastDetSvc", detSvc);
  if (sc.isFailure() ) {
    msglog << MSG::ERROR << "  can't get GlastDetSvc " << endreq;
    return sc;
  }

  //-- Get GlastDetSvc constants --//
  PARAMAP param;
  param[&m_maxAdc]       = string("cal.maxAdcValue");

  for(PARAMAP::iterator it=param.begin(); it!=param.end();it++)
    if(!detSvc->getNumericConstByName((*it).second, &val)) {
      msglog << MSG::ERROR << " constant " <<(*it).second <<" not defiPed" << endreq;
      return StatusCode::FAILURE;
    } else *((*it).first)=(int)val;

  
  // obtain CalCalibSvc
  sc = service(m_calCalibSvcName.value(), m_calCalibSvc);
  if (sc.isFailure()) {
    msglog << MSG::ERROR << "can't get " << m_calCalibSvcName  << endreq;
    return sc;
  }

  // open optional tuple file
  if (m_tupleFilename.value().length() > 0 ) {
    m_tupleFile.reset(new TFile(m_tupleFilename.value().c_str(),"RECREATE","XtalDigiTuple"));
    if (!m_tupleFile.get())
      // allow to continue w/out tuple file as it is not a _real_ failure
      msglog << MSG::ERROR << "Unable to create TTree object: " << m_tupleFilename << endreq;
    
    else {
      m_tuple = new TTree("XtalDigiTuple","XtalDigiTuple");
      if (!m_tuple) {
        msglog << MSG::ERROR << "Unable to create tuple" << endreq;
        return StatusCode::FAILURE;
      }

      //-- Add Branches to tree --//
      if (!m_tuple->Branch("RunID",          &m_dat.RunID,             "RunID/i")            ||
          !m_tuple->Branch("EventID",        &m_dat.EventID,           "EventID/i")          ||
          !m_tuple->Branch("XtalIdx",        &m_dat.xtalIdx,            "xtalIdx/i")          ||
          !m_tuple->Branch("adcPed",         m_dat.adcPed.begin(),     "adcPed[2][4]/F")     ||
          !m_tuple->Branch("diodeCIDAC",      m_dat.diodeCIDAC.begin(),"diodeCIDAC[2][2]/F") ||
          !m_tuple->Branch("ped",            m_dat.ped.begin(),        "ped[2][4]/F")        ||
          !m_tuple->Branch("lacThreshCIDAC", m_dat.lacThreshCIDAC.begin(),  "lacThreshCIDAC[2]/F")     ||
          !m_tuple->Branch("uldTholdADC",    m_dat.uldTholdADC.begin(),   "uldTholdADC[2][4]/F")   ||
          !m_tuple->Branch("rng",            m_dat.rng.begin(),        "rng[2]/b")           ||
          !m_tuple->Branch("lac",            m_dat.lac.begin(),        "lac[2]/b")           ||
          !m_tuple->Branch("saturated",      m_dat.saturated.begin(),  "saturated[2]/b")) {
        msglog << MSG::ERROR << "Couldn't create tuple branch" << endreq;
        return StatusCode::FAILURE;
      }
    }
 
    //-- Retreive EventDataSvc (used by tuple)
    sc = serviceLocator()->service( "EventDataSvc", m_evtSvc, true );
    if(sc.isFailure()){
      msglog << MSG::ERROR << "Could not find EventDataSvc" << endreq;
      return sc;
    }


  } // optional tuple

  ///////////////////////////////////////
  //-- RETRIEVE HELPER TOOLS & SVCS  --//
  ///////////////////////////////////////


  //-- find optional CalFailureModeSvc --//
  sc = service("CalFailureModeSvc", m_calFailureModeSvc);
  if (sc.isFailure() ) {
    msglog << MSG::INFO << "  Did not find CalFailureMode service" << endreq;
    m_calFailureModeSvc = 0;
  }

  // this tool may also be shared by CalTrigTool, global ownership
  sc = toolSvc()->retrieveTool("PrecalcCalibTool", 
                               m_precalcCalibName, 
                               m_precalcCalib,
                               0); // shared by other code
  if (sc.isFailure() ) {
    msglog << MSG::ERROR << "  Unable to create " << m_precalcCalibName << endreq;
    return sc;
  }

  return StatusCode::SUCCESS;
}

/** \brief calculate single cystal CalDigi object from diode signal levels

  Basic Algorithm:
  - test LAC thresholds, quit early if zeroSuppression is enabled
  - check channel failure status from CalFailureMode
  - convert cidac->adc units for each channel
  - select 'best range' or lowest non-saturated range
  - generate CalDigi object
  - optionally fill debugging tuple.
  
*/
StatusCode XtalDigiTool::calculate(const ICalSignalTool::XtalSignalMap &cidac,
                                   Event::CalDigi &calDigi,
                                   CalUtil::CalArray<CalUtil::FaceNum, bool> &lacBits,
                                   const bool zeroSuppress) {
  StatusCode sc;

  m_dat.Clear();

  const CalXtalId xtalId = calDigi.getPackedId();
  m_dat.xtalIdx = XtalIdx(xtalId);
  
  ////////////////////////
  // Stage 1: LAC TESTS //
  ////////////////////////
    
  for (FaceNum face; face.isValid(); face++) {
    const FaceIdx faceIdx(m_dat.xtalIdx, face);
    sc = m_precalcCalib->getLacCIDAC(faceIdx, m_dat.lacThreshCIDAC[face]);
    if (sc.isFailure()) return sc;
    
    // set log-accept flags
    lacBits[face] = cidac[XtalDiode(face, LRG_DIODE)] >= 
      m_dat.lacThreshCIDAC[face];
  }

  //-- Quick Exit --//
  // once LAC & triggers have been processed (in CIDAC scale
  // I can now exit and save rest of processing including adc->dac
  // calculations
  if (zeroSuppress && !lacBits[POS_FACE] && !lacBits[NEG_FACE])
    return StatusCode::SUCCESS;

  /////////////////////////////////////////
  // Stage 2: Optional failure mode test //
  /////////////////////////////////////////
  sc = getFailureStatus(lacBits);
  if (sc.isFailure()) return sc;
  
  ///////////////////////////////
  // Stage 3: convert cidac->adc //
  ///////////////////////////////
  for (XtalRng xRng; xRng.isValid(); xRng++) {
    const XtalDiode xDiode(xRng.getFace(), xRng.getRng().getDiode());
    // calcuate adc val from cidac
    sc = m_calCalibSvc->evalADC(RngIdx(m_dat.xtalIdx, xRng), 
                                cidac[xDiode],
                                m_dat.adcPed[xRng]);
    if (sc.isFailure()) return sc;   
  } // xRng
  
  //////////////////////////////
  // Stage 4: Range Selection //
  //////////////////////////////
  sc = rangeSelect();
  if (sc.isFailure()) return sc;
    
  ///////////////////////////////////////
  //-- STEP 5: Populate Return Vars --//
  ///////////////////////////////////////
  // generate xtalDigReadouts.
  sc = fillDigi(calDigi);
  if (sc.isFailure()) return sc;

  /////////////////////////////////////////////////////////
  //-- STEP 6: Populate XtalDigiTuple vars (optional) --//
  /////////////////////////////////////////////////////////
  
  // following steps only needed if tuple output is selected
  if (m_tuple) {
    // if lac_only, then only continue if one of 2 lac flags is high
    if ((m_tupleLACOnly && (lacBits[NEG_FACE] || 
                            lacBits[POS_FACE])) ||
        !m_tupleLACOnly) {

      // get pointer to EventHeader (has runId, evtId, etc...)
      Event::EventHeader *evtHdr = 0;
      if (m_evtSvc) {
        evtHdr = SmartDataPtr<Event::EventHeader>(m_evtSvc, EventModel::EventHeader) ;
        if (!evtHdr) {
          MsgStream msglog( msgSvc(), name() );
          msglog << MSG::ERROR << "Event header not found !" << endreq ;
        }

        m_dat.RunID   = evtHdr->run();
        m_dat.EventID = evtHdr->event();
      }
      
      // populate cidac fields
      copy(cidac.begin(),
           cidac.end(),
           m_dat.diodeCIDAC.begin());
      
      // populate lac bits
      // wack MSVC warning doesn't like me cast bool to int.
      m_dat.lac[POS_FACE] = lacBits[POS_FACE] != 0;
      m_dat.lac[NEG_FACE] = lacBits[NEG_FACE] != 0;

      // instruct tuple to fill
      m_tuple->Fill();
    }
    
  }

  return StatusCode::SUCCESS;
}

/// select lowest non-saturated ADC range on each crystal face
StatusCode XtalDigiTool::rangeSelect() {
  for (FaceNum face; face.isValid(); face++) {
    // retrieve calibration

    //-- TholdCI --//
    const FaceIdx faceIdx(m_dat.xtalIdx, face);
    CalTholdCI const*const tholdCI = m_calCalibSvc->getTholdCI(faceIdx);
    if (!tholdCI) return StatusCode::FAILURE;;

    // BEST RANGE
    RngNum rng;
    for (rng=LEX8; rng.isValid(); rng++) {
      // get ULD threshold
      XtalRng xRng(face,rng);
      m_dat.uldTholdADC[xRng] = tholdCI->getULD(rng.val())->getVal();

      // case of HEX1 range 
      // adc vals have a ceiling in HEX1
      if (rng == HEX1) {
        if (m_dat.adcPed[xRng] >= m_dat.uldTholdADC[xRng]) {
          // set ADC to max val
          m_dat.adcPed[xRng] =  m_dat.uldTholdADC[xRng];
          
          // set 'pegged' flag
          m_dat.saturated[face] = true;
        }
        
        // break before rng is incremented out-of-bounds
        break; 
      } else { // 1st 3 ranges
        // break on 1st energy rng that is < threshold
        if (m_dat.adcPed[xRng] <= m_dat.uldTholdADC[xRng]) break;
      }
    }
    
    // assign range selection
    m_dat.rng[face] = rng;
  }  // per face, range selection

  return StatusCode::SUCCESS;
}

// fill digi with adc values
StatusCode XtalDigiTool::fillDigi(CalDigi &calDigi) {
  
  //-- How many readouts ? --//
  int roLimit;
  const short rangeMode = calDigi.getMode();
  switch (rangeMode) {
  case CalXtalId::BESTRANGE:
    roLimit = 1;
    break;
  case CalXtalId::ALLRANGE:
    roLimit = 4;
    break;
  default:
    MsgStream msglog(msgSvc(), name()); 
    msglog << MSG::ERROR << "Unsupported Cal readout Mode" << endreq;
    return StatusCode::FAILURE;
  }

  // set up the digi
  CalArray<FaceNum, RngNum> roRange;
  CalArray<FaceNum, float> adc;
  for (int nRo=0; nRo < roLimit; nRo++) {
    for (FaceNum face; face.isValid(); face++) {
      // represents ranges used for current readout in loop
      roRange[face] = RngNum((m_dat.rng[face].val() + nRo) % RngNum::N_VALS); 

      const XtalRng xRng(face,roRange[face]);
      const RngIdx rngIdx(m_dat.xtalIdx, face, roRange[face]);
      CalibData::Ped const*const ped = m_calCalibSvc->getPed(rngIdx);
      if (!ped) return StatusCode::FAILURE;
      m_dat.ped[xRng] = ped->getAvr();
  
      ////////////////////////////////////////
      // Stage 5: ADD PEDS, CHECK ADC RANGE //
      ////////////////////////////////////////
      // double check that ADC vals are all >= 0
      // double check that ADC vals are all < maxadc(4095)
  
      // must be after rangeSelect()  b/c rangeSelect()
      // may clip HEX1 to HEX1 saturation point
      adc[face] = max<float>(0, m_dat.adcPed[xRng] + m_dat.ped[xRng]);
      adc[face] = round_int(min<float>(m_maxAdc, adc[face]));
    }
      
    CalDigi::CalXtalReadout ro = CalDigi::CalXtalReadout(roRange[POS_FACE].val(), 
                                                         (short)adc[POS_FACE], 
                                                         roRange[NEG_FACE].val(), 
                                                         (short)adc[NEG_FACE], 
                                                         m_dat.failureStatus);
    calDigi.addReadout(ro);
  }
  
  
  return StatusCode::SUCCESS;
}

StatusCode XtalDigiTool::finalize() {
  // make sure optional tuple is closed out                                                        
  if (m_tupleFile.get()) {
    m_tupleFile->Write();
    m_tupleFile->Close(); // trees deleted                                                         
  }

  return StatusCode::SUCCESS;
}
  
StatusCode XtalDigiTool::getFailureStatus(const CalArray<FaceNum, bool> &lacBits) {
  /// quick exit if option not enabled
  if (m_calFailureModeSvc == 0)
    return StatusCode::SUCCESS;
  
  const CalXtalId xtalId(m_dat.xtalIdx.getCalXtalId());

  if (m_calFailureModeSvc->matchChannel(xtalId,
                                        (CalXtalId::POS))) {
    
    if (lacBits[POS_FACE]) (m_dat.failureStatus |= Event::CalDigi::CalXtalReadout::DEAD_P);
  }
  if (m_calFailureModeSvc->matchChannel(xtalId,
                                        (CalXtalId::NEG))) {
    if (lacBits[NEG_FACE]) (m_dat.failureStatus |= Event::CalDigi::CalXtalReadout::DEAD_N);
    
  }

  if ((m_dat.failureStatus & 0x00FF) == 0) m_dat.failureStatus |= Event::CalDigi::CalXtalReadout::OK_P;
  if ((m_dat.failureStatus & 0xFF00) == 0) m_dat.failureStatus |= Event::CalDigi::CalXtalReadout::OK_N;

  return StatusCode::SUCCESS;
}


