// LOCAL
#include "CalXtalResponse/IXtalEneTool.h"
#include "CalXtalResponse/ICalCalibSvc.h"

// GLAST
#include "CalUtil/CalDefs.h"

// EXTLIB
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/MsgStream.h"

// STD
#include <cmath>

using namespace CalDefs;

/*! @class XtalEneTool
  \author Zachary Fewtrell
  \brief Simple implementation of IXtalEneTool.  Faithfully pasted from XtalRecAlg v5r6p1
*/

class XtalEneTool : 
  public AlgTool, 
  virtual public IXtalEneTool 
{
public:

  /// default ctor, declares jobOptions.
  XtalEneTool::XtalEneTool( const string& type, 
                            const string& name, 
                            const IInterface* parent);

  /// retrieves needed parameters and pointers to required services
  virtual StatusCode initialize();

  /// calculate energy deposition given the digi response for both xtal faces
  StatusCode calculate(const CalXtalId &xtalId, 
                       CalXtalId::AdcRange rngP,
                       CalXtalId::AdcRange rngN,
                       int adcP, 
                       int adcN,
                       float &energy,     
                       bool  &belowThresh 
                       );

  /// calculate energy deposition given the digi response for one xtal face
  StatusCode calculate(const CalXtalId &xtalId,
                       int adc, 
                       float position,
                       float &energy,     
                       bool  &belowThresh 
                       );
private:

  StringProperty m_calCalibSvcName;                         ///< name of CalCalibSvc to use for calib constants.
  ICalCalibSvc *m_calCalibSvc;                             ///< pointer to CalCalibSvc object.

};

static ToolFactory<XtalEneTool> s_factory;
const IToolFactory& XtalEneToolFactory = s_factory;

XtalEneTool::XtalEneTool( const string& type,
                          const string& name,
                          const IInterface* parent)
  : AlgTool(type,name,parent) {
  declareInterface<IXtalEneTool>(this);

  declareProperty("CalCalibSvc", m_calCalibSvcName="CalCalibSvc");
}

StatusCode XtalEneTool::initialize() {
  MsgStream msglog(msgSvc(), name());
  msglog << MSG::INFO << "initialize" << endreq;

  StatusCode sc;

  // obtain CalCalibSvc
  sc = service(m_calCalibSvcName.value(), m_calCalibSvc);
  if (sc.isFailure()) {
    msglog << MSG::ERROR << "Unable to get CalCalibSvc." << endreq;
    return sc;
  }

  return StatusCode::SUCCESS;
}

StatusCode XtalEneTool::calculate(const CalXtalId &xtalId, 
                                  CalXtalId::AdcRange rngP,
                                  CalXtalId::AdcRange rngN,
                                  int adcP, 
                                  int adcN,
                                  float &energy,
                                  bool &belowThresh
                                  ) {
  StatusCode sc;
  XtalIdx xtalIdx(xtalId); // used for array indexing

  // initial return vals
  energy = 0;
  belowThresh = false;
 

  //-- RETRIEVE PEDESTAL (POS_FACE) --//
  float pedP, pedN, sig, cos;
  RngIdx rngIdxP(xtalIdx, POS_FACE, rngP);
  sc = m_calCalibSvc->getPed(rngIdxP.getCalXtalId(), pedP, sig, cos);
  if (sc.isFailure()) return sc;
  int adcPedP = adcP - pedP;   // ped subtracted ADC
  
  //-- THROW OUT LOW ADC VALS (POS_FACE) --/
  CalibData::ValSig fle,fhe,lacP, lacN;
  CalXtalId tmpIdP(xtalId.getTower(),
                   xtalId.getLayer(),
                   xtalId.getColumn(),
                   POS_FACE);
  sc = m_calCalibSvc->getTholdCI(tmpIdP,fle,fhe,lacP);
  if (adcPedP < lacP.getVal()*0.5) {
    belowThresh = true;
    return StatusCode::SUCCESS;
  }

  //-- RETRIEVE PEDESTAL (NEG_FACE) --//
  RngIdx rngIdxN(xtalIdx, NEG_FACE, rngN);
  sc = m_calCalibSvc->getPed(rngIdxN.getCalXtalId(), pedN, sig, cos);
  if (sc.isFailure()) return sc;
  int adcPedN = adcN - pedN;

  //-- THROW OUT LOW ADC VALS (NEG_FACE) --/
  CalXtalId tmpIdN(xtalId.getTower(),
                   xtalId.getLayer(),
                   xtalId.getColumn(),
                   NEG_FACE);
  sc = m_calCalibSvc->getTholdCI(tmpIdN,fle,fhe,lacN);
  if (adcPedN < lacN.getVal()*0.5) {
    belowThresh = true;
    return StatusCode::SUCCESS;
  }
  

  //-- CONVERT ADCs -> DAC --//
  double dacP, dacN;
  sc = m_calCalibSvc->evalDAC(rngIdxP.getCalXtalId(), adcPedP, dacP);
  if (sc.isFailure()) return sc;
  sc = m_calCalibSvc->evalDAC(rngIdxN.getCalXtalId(), adcPedN, dacN);
  if (sc.isFailure()) return sc;

  
  //-- CONVERT DIODE SIZE (IF NEEDED) --//
  // if diodes are different, then I need to convert them both to sm diode
  DiodeNum diodeP = RngNum(rngP).getDiode();
  DiodeNum diodeN = RngNum(rngN).getDiode();
  if (diodeP != diodeN) { 
    // get small diode asymetry at center of xtal (used for both conversions)
    // NOTE:
    // using center of xtal for position since the lrg/sm ratio is 
    // effectively constant throughout the xtal length and this keeps the
    // energy calculation independent of the position calculation.
    double asymSm;
    sc = m_calCalibSvc->evalAsymSm(xtalId, 0.0, asymSm);
    if (sc.isFailure()) return sc;

    if (diodeP == LRG_DIODE) {
      // STRATEGY:
      // convert PL 2 PS
      // AsymNSPB = log(PL / NS)
      // AsymSm = log(PS / NS)
      // exp(sm)/exp(nspb) = (PS/NS)/(PL/NS) = PS/PL
      // exp(nspb-sm) = PS/PL
      
      double asymNSPB;
      sc = m_calCalibSvc->evalAsymNSPB(xtalId, 0.0, asymNSPB);
      if (sc.isFailure()) return sc;

      dacP /= exp(asymSm - asymNSPB);
      diodeP = SM_DIODE;
    } else { // case: diodeN == LRG
      // STRATEGY:
      // convert NL -> NS
      // asymPSNB = log(PS/NL)
      // asymSm = log(PS/NS)
      // exp(psnb)/exp(sm) = (PS/NL)/(PS/NS)
      // exp(psnb-sm) = NS/NL

      double asymPSNB;
      sc = m_calCalibSvc->evalAsymPSNB(xtalId, 0.0, asymPSNB);
      if (sc.isFailure()) return sc;

      dacN *= exp(asymPSNB-asymSm);
      diodeN = SM_DIODE;
    }
  }
  // now both dac values are for the same diode size

  //-- RETRIEVE MEVPERDAC --//
  CalibData::ValSig mpdLrg, mpdSm;
  sc = m_calCalibSvc->getMeVPerDac(xtalId, mpdLrg, mpdSm);
  if (sc.isFailure()) return sc;
  
  //-- CALC MEAN DAC & ENERGY --//
  double meanDAC = sqrt(dacP*dacN);
  if (diodeP == SM_DIODE)
    energy = meanDAC*mpdSm.getVal();
  else
    energy = meanDAC*mpdLrg.getVal();

  return StatusCode::SUCCESS;
}

StatusCode XtalEneTool::calculate(const CalXtalId &xtalId,
                                  int adc, 
                                  float position,
                                  float &energy,  
                                  bool &belowThresh
                                  ) {
  StatusCode sc;

  energy = 0;
  belowThresh = false;

  // check that CalXtalId has face & range fields enabled
  if (!xtalId.validFace() || !xtalId.validRange())
    throw invalid_argument("XtalEneTool: CalXtalId requires valid"
                           " face & range fields");

  // retrieve pedestal
  float ped, sig, cos;
  sc = m_calCalibSvc->getPed(xtalId, ped, sig, cos);
  if (sc.isFailure()) return sc;

  //-- THROW OUT LOW ADCS --/
  CalibData::ValSig fle,fhe,lac;

  // retreive threshold
  sc = m_calCalibSvc->getTholdCI(xtalId,fle,fhe,lac);

  // set flag
  if (adc < lac.getVal()*0.5) {
    belowThresh = true;
    return StatusCode::SUCCESS;
  }

  // convert adc->dac units
  double dac;
  sc = m_calCalibSvc->evalDAC(xtalId, adc-ped, dac);
  if (sc.isFailure()) return sc;

  RngNum rng(xtalId.getRange());
  DiodeNum diode(rng.getDiode());

  // use appropriate asymmetry to calc position based on
  // diode size
  double asym;
  if (diode == SM_DIODE) {
    sc = m_calCalibSvc->evalAsymSm(xtalId, position, asym);
    if (sc.isFailure()) return sc;
  } else {
    sc = m_calCalibSvc->evalAsymLrg(xtalId, position, asym);
    if (sc.isFailure()) return sc;
  }

  // calculate the 'other' dac
  FaceNum face = xtalId.getFace();
  double otherDac;
  if (face == POS_FACE) otherDac = dac/exp(asym);
  else otherDac = dac*exp(asym);

  double meanDAC = sqrt(dac*otherDac);

  CalibData::ValSig mpdLrg, mpdSm;
  sc = m_calCalibSvc->getMeVPerDac(xtalId, mpdLrg, mpdSm);
  if (sc.isFailure()) return sc;

  if (diode == SM_DIODE)
    energy = meanDAC*mpdSm.getVal();
  else
    energy = meanDAC*mpdLrg.getVal();

  
  return StatusCode::SUCCESS;
}
