// Include files

// LOCAL
#include "CalXtalResponse/IXtalPosTool.h"
#include "CalXtalResponse/ICalCalibSvc.h"

// GLAST
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "CalUtil/CalDefs.h"

// EXTLIB
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/MsgStream.h"

// STD
#include <cmath>

using namespace CalDefs;

/*! @class XtalPosTool
  \author Zachary Fewtrell
  \brief ICalPosTool implementation uses calibGenCAL v3.0 constants.
*/
class XtalPosTool : public AlgTool, virtual public IXtalPosTool {
public:

  /// default ctor, declares jobOptions
  XtalPosTool::XtalPosTool( const string& type, 
                            const string& name, 
                            const IInterface* parent);

  /// retrieves needed paramters and pointers to required services
  virtual StatusCode initialize();

  /// calculate position of energy deposition in mm from xtal-center = 0
  StatusCode calculate(const CalXtalId &xtalId,
                       CalXtalId::AdcRange rngP,
                       CalXtalId::AdcRange rngN,
                       int adcP, 
                       int adcN, 
                       float &position                // output
                       );
private:

  StringProperty m_calCalibSvcName;                   ///< name of CalCalibSvc to use for calib constants.
  ICalCalibSvc *m_calCalibSvc;                        ///< pointer to CalCalibSvc object.

  double m_CsILength;                    ///< Xtal length
};

static ToolFactory<XtalPosTool> s_factory;
const IToolFactory& XtalPosToolFactory = s_factory;

XtalPosTool::XtalPosTool( const string& type,
                          const string& name,
                          const IInterface* parent)
  : AlgTool(type,name,parent) {
  declareInterface<IXtalPosTool>(this);

  declareProperty("CalCalibSvc", m_calCalibSvcName="CalCalibSvc");
}

StatusCode XtalPosTool::initialize() {
  MsgStream msglog(msgSvc(), name());
  msglog << MSG::INFO << "initialize" << endreq;

  StatusCode sc;

  // obtain CalCalibSvc
  sc = service(m_calCalibSvcName.value(), m_calCalibSvc);
  if (sc.isFailure()) {
    msglog << MSG::ERROR << "Unable to get CalCalibSvc." << endreq;
    return sc;
  }

  // try to find the GlastDevSvc service
  IGlastDetSvc* detSvc;
  sc = service("GlastDetSvc", detSvc);
  if (sc.isFailure() ) {
    msglog << MSG::ERROR << "  Unable to retrieve GlastDetSvc " << endreq;
    return sc;
  }

  // doubles are done separately
  sc = detSvc->getNumericConstByName("CsILength", &m_CsILength);
  if (sc.isFailure()) {
    msglog << MSG::ERROR << " constant CsILength not defined" << endreq;
    return StatusCode::FAILURE;
  }

  return StatusCode::SUCCESS;
}

StatusCode XtalPosTool::calculate(const CalXtalId &xtalId,
                                  CalXtalId::AdcRange rngP,
                                  CalXtalId::AdcRange rngN,
                                  int adcP, 
                                  int adcN, 
                                  float &position                // output
                                  ) {
  StatusCode sc;

  position = 0;

  // used for index manipulation
  XtalIdx xtalIdx(xtalId);

  //-- RETRIEVE PED (POS_FACE) --//
  float pedP, pedN, sig, cos;
  RngIdx rngIdxP(xtalIdx, POS_FACE, rngP);
  sc = m_calCalibSvc->getPed(rngIdxP.getCalXtalId(), pedP, sig, cos);
  if (sc.isFailure()) return sc;
  float adcPedP = adcP - pedP;   // ped subtracted ADC
  
  //-- RETRIEVE PED (NEG_FACE) --//
  RngIdx rngIdxN(xtalIdx, NEG_FACE, rngN);
  sc = m_calCalibSvc->getPed(rngIdxN.getCalXtalId(), pedN, sig, cos);
  if (sc.isFailure()) return sc;
  float adcPedN = adcN - pedN;
    
  // convert adc->dac units
  double dacP, dacN;
  sc = m_calCalibSvc->evalDAC(rngIdxP.getCalXtalId(), adcPedP, dacP);
  if (sc.isFailure()) return sc;
  sc = m_calCalibSvc->evalDAC(rngIdxN.getCalXtalId(), adcPedN, dacN);
  if (sc.isFailure()) return sc;

  // check for invalid dac vals (i need to take the logarithm)
  if (dacP <= 0 || dacN <=0) {
    // create MsgStream only when needed for performance
    MsgStream msglog(msgSvc(), name()); 
    msglog << MSG::WARNING << "DAC val <= 0, can't calculate position.  This shouldn't happen." << endl;
    return StatusCode::FAILURE;
  }
  double asym = log(dacP/dacN);

  // select proper asymmetry function & find position
  double pos;
  DiodeNum diodeP = RngNum(rngP).getDiode();
  DiodeNum diodeN = RngNum(rngN).getDiode();
  if (diodeP == LRG_DIODE) {
    if (diodeN == LRG_DIODE)
      sc = m_calCalibSvc->evalPosLrg(xtalId, asym, pos);
    else
      sc = m_calCalibSvc->evalPosNSPB(xtalId, asym, pos);
  } else {
    if (diodeN == SM_DIODE)
      sc = m_calCalibSvc->evalPosSm(xtalId, asym, pos);
    else
      sc = m_calCalibSvc->evalPosPSNB(xtalId, asym, pos);
  }
  if (sc.isFailure()) return sc;

  // report any pos which evals to outside the xtal
  // as on the edge of the xtal
  
  double halfXtalLen = m_CsILength/2;
  pos = min(halfXtalLen, pos);
  pos = max(-1*halfXtalLen, pos);

  position = pos;

  return StatusCode::SUCCESS;
}
