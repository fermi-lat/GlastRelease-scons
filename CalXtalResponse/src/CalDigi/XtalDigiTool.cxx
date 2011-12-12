//    $Header$

/** @file     
    @author Z.Fewtrell

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
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/IDataProviderSvc.h"

// STD
#include <map>
#include <cmath>

using namespace CalUtil;
using namespace Event;
using namespace CalibData;
using namespace idents;
using namespace std;

//static ToolFactory<XtalDigiTool> s_factory;
//const IToolFactory& XtalDigiToolFactory = s_factory;
DECLARE_TOOL_FACTORY(XtalDigiTool);

static float round_int(float in) { return floor(in + 0.5);}

XtalDigiTool::XtalDigiTool( const string& type,
                            const string& name,
                            const IInterface* parent)
  : AlgTool(type,name,parent),
    m_calCalibSvc(0),
    m_maxAdc(-1),
    m_calFailureModeSvc(0),
    m_precalcCalib(0),
    m_calSignalTool(0)
{
  declareInterface<IXtalDigiTool>(this);

  declareProperty("CalCalibSvc",        m_calCalibSvcName    = "CalCalibSvc");
  declareProperty("PrecalcCalibTool",   m_precalcCalibName   = "PrecalcCalibTool");
  declareProperty("CalSignalToolName",   m_calSignalToolName = "CalSignalTool");


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

  sc = toolSvc()->retrieveTool("CalSignalTool",
                               m_calSignalToolName,
                               m_calSignalTool,
                               0); // intended to be shared
  if (sc.isFailure() ) {
    msglog << MSG::ERROR << "  can't create " << m_calSignalToolName << endreq;
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
  
*/

StatusCode XtalDigiTool::calculate(Event::CalDigi &calDigi,
                                   CalVec<FaceNum, bool> &lacBits,
                                   const bool zeroSuppress, string calFirstRng) {
  StatusCode sc;
  MsgStream msglog(msgSvc(), name());
  
  const CalXtalId xtalId = calDigi.getPackedId();
  const XtalIdx xtalIdx(xtalId);

  ////////////////////////////////////////////////////////
  //-- Stage 0: retrieve signal levels for each diode --//
  ////////////////////////////////////////////////////////
  CalArray<XtalDiode, float> cidac;
  for (XtalDiode xDiode; xDiode.isValid(); xDiode++) {
    sc = m_calSignalTool->getDiodeSignal(DiodeIdx(xtalIdx, xDiode), cidac[xDiode]);
    if (sc.isFailure())
      return sc;
  }
  
  ////////////////////////
  // Stage 1: LAC TESTS //
  ////////////////////////
  
  for (FaceNum face; face.isValid(); face++) {
    const FaceIdx faceIdx(xtalIdx, face);
    float lacThreshCIDAC;
    sc = m_precalcCalib->getLacCIDAC(faceIdx, lacThreshCIDAC);
    if (sc.isFailure()) return sc;
    
    // set log-accept flags
    lacBits[face] = cidac[XtalDiode(face, LRG_DIODE)] >= 
      lacThreshCIDAC;
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
  unsigned short failureStatus;
  sc = getFailureStatus(xtalId, lacBits, failureStatus);
  if (sc.isFailure()) return sc;
  
  ///////////////////////////////
  // Stage 3: convert cidac->adc //
  ///////////////////////////////
  CalArray<XtalRng, float> adcPed;
  for (XtalRng xRng; xRng.isValid(); xRng++) {
    const XtalDiode xDiode(xRng.getFace(), xRng.getRng().getDiode());
    // calcuate adc val from cidac
    const RngIdx rngIdx(xtalIdx, xRng);
    sc = m_calCalibSvc->evalADC(rngIdx, 
                                cidac[xDiode],
                                adcPed[xRng]);
    if (sc.isFailure()) return sc;   

#if 0
    cout << "XtalDigiTool: " << rngIdx.toStr()
         << " cidac " << cidac[xDiode]
         << " adcPed " << adcPed[xRng]
         << endl;
#endif

  } // xRng
  
  //////////////////////////////
  // Stage 4: Range Selection //
  //////////////////////////////

  ///////////////////////////////////////
  //-- STEP 5: Populate Return Vars --//
  ///////////////////////////////////////
  
  // generate xtalDigReadouts.

  if(calFirstRng== "autoRng")  // default (best range) R/O mode
    {
      CalVec<FaceNum, RngNum> bestRng;
      
      sc= rangeSelect(xtalIdx, adcPed, bestRng);
      if(sc.isFailure()) return sc;

      sc= fillDigi(calDigi, adcPed, bestRng, failureStatus);
    }
  else                         // forced range R/O mode
    {
      if(calFirstRng== "lex8" || calFirstRng== "lex1" || calFirstRng== "hex8" || calFirstRng== "hex1")
        {
          CalVec<FaceNum, RngNum> forcRng;
          
          if(calFirstRng== "lex8")
            {
              forcRng[POS_FACE]= LEX8;
              forcRng[NEG_FACE]= LEX8;
            }
          
          if(calFirstRng== "lex1")
            {
              forcRng[POS_FACE]= LEX1;
              forcRng[NEG_FACE]= LEX1;
            }
          
          if(calFirstRng== "hex8")
            {
              forcRng[POS_FACE]= HEX8;
              forcRng[NEG_FACE]= HEX8;
            }
          
          if(calFirstRng== "hex1")
            {
              forcRng[POS_FACE]= HEX1;
              forcRng[NEG_FACE]= HEX1;
            }
          
          sc= fillDigi(calDigi, adcPed, forcRng, failureStatus);
        }
      else
        {
          msglog << MSG::ERROR << "invalid R/O range specified in jobOptions" << endreq;
          return StatusCode::FAILURE;
        }
    }
      
  if(sc.isFailure()) return sc;

  return StatusCode::SUCCESS;
}

/// select lowest non-saturated ADC range on each crystal face
StatusCode XtalDigiTool::rangeSelect(const XtalIdx xtalIdx,
                                     CalArray<XtalRng, float> &adcPed,
                                     CalVec<FaceNum, RngNum> &bestRng) {
  for (FaceNum face; face.isValid(); face++) {
    // retrieve calibration

    //-- TholdCI --//
    const FaceIdx faceIdx(xtalIdx, face);
    CalTholdCI const*const tholdCI = m_calCalibSvc->getTholdCI(faceIdx);
    if (!tholdCI) return StatusCode::FAILURE;;

    // BEST RANGE
    RngNum rng;
    for (rng=LEX8; rng.isValid(); rng++) {
      // get ULD threshold
      const XtalRng xRng(face,rng);
      const float uldTholdADC = tholdCI->getULD(rng.val())->getVal();

      // case of HEX1 range 
      // adc vals have a ceiling in HEX1
      if (rng == HEX1) {
        if (adcPed[xRng] >= uldTholdADC)
          adcPed[xRng] =  uldTholdADC;      // set ADC to max val
        
        // break before rng is incremented out-of-bounds
        break; 

      } else  // non-HEX1 case
        // break on 1st energy rng that is < threshold
        if (adcPed[xRng] <= uldTholdADC) break;

    }
    
    // assign range selection
    bestRng[face] = rng;
  }  // per face, range selection

  return StatusCode::SUCCESS;
}

// fill digi with adc values
StatusCode XtalDigiTool::fillDigi(CalDigi &calDigi,
                                  const CalArray<XtalRng, float> &adcPed,
                                  const CalVec<FaceNum, RngNum> &bestRng,
                                  unsigned short failureStatus
                                  ) {
  
  //-- How many readouts ? --//
  unsigned short roLimit;
  const unsigned short rangeMode = calDigi.getMode();
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
  CalVec<FaceNum, RngNum> roRange;
  CalVec<FaceNum, float> adc;
  const CalXtalId xtalId(calDigi.getPackedId());
  for (unsigned short nRo=0; nRo < roLimit; nRo++) {
    for (FaceNum face; face.isValid(); face++) {
      // represents ranges used for current readout in loop
      roRange[face] = RngNum((bestRng[face].val() + nRo) % RngNum::N_VALS); 

      const XtalRng xRng(face,roRange[face]);
      const XtalIdx xtalIdx(xtalId);
      const RngIdx rngIdx(xtalIdx, face, roRange[face]);
      float ped;
      StatusCode sc = m_calCalibSvc->getPed(rngIdx,ped);
      if (sc.isFailure()) return StatusCode::FAILURE;

  
      ////////////////////////////////////////
      // Stage 5: ADD PEDS, CHECK ADC RANGE //
      ////////////////////////////////////////
      // double check that ADC vals are all >= 0
      // double check that ADC vals are all < maxadc(4095)
  
      // must be after rangeSelect()  b/c rangeSelect()
      // may clip HEX1 to HEX1 saturation point
      adc[face] = max<float>(0, adcPed[xRng] + ped);
      adc[face] = round_int(min<float>(m_maxAdc, adc[face]));

#if 0
      cout << "fillDigi " << rngIdx.toStr()
           << " adcPed " << adcPed[xRng]
           << " ped " << ped
           << " adc " << adc[face]
           << endl;
#endif

    }
      
    CalDigi::CalXtalReadout ro = CalDigi::CalXtalReadout(roRange[POS_FACE].val(), 
                                                         (short)adc[POS_FACE], 
                                                         roRange[NEG_FACE].val(), 
                                                         (short)adc[NEG_FACE], 
                                                         failureStatus);

    calDigi.addReadout(ro);
  }
  
  
  return StatusCode::SUCCESS;
}

StatusCode XtalDigiTool::finalize() {
  return StatusCode::SUCCESS;
}
  
StatusCode XtalDigiTool::getFailureStatus(const idents::CalXtalId xtalId,
                                          const CalVec<FaceNum, bool> &lacBits,
                                          unsigned short &failureStatus
                                          ) {
  /// quick exit if option not enabled
  if (m_calFailureModeSvc == 0)
    return StatusCode::SUCCESS;
  
  if (m_calFailureModeSvc->matchChannel(xtalId,
                                        (CalXtalId::POS))) {
    
    if (lacBits[POS_FACE]) (failureStatus |= Event::CalDigi::CalXtalReadout::DEAD_P);
  }
  if (m_calFailureModeSvc->matchChannel(xtalId,
                                        (CalXtalId::NEG))) {
    if (lacBits[NEG_FACE]) (failureStatus |= Event::CalDigi::CalXtalReadout::DEAD_N);
    
  }

  if ((failureStatus & 0x00FF) == 0) failureStatus |= Event::CalDigi::CalXtalReadout::OK_P;
  if ((failureStatus & 0xFF00) == 0) failureStatus |= Event::CalDigi::CalXtalReadout::OK_N;

  return StatusCode::SUCCESS;
}


