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

  energy = 0;

  // used for index manipulation
  XtalIdx xtalIdx(xtalId);

  // retrieve pedestals
  RngIdx rngIdxP(xtalIdx, POS_FACE, rngP);
  RngIdx rngIdxN(xtalIdx, NEG_FACE, rngN);
  
  float pedP, pedN, sig, cos;
  sc = m_calCalibSvc->getPed(rngIdxP.getCalXtalId(), pedP, sig, cos);
  if (sc.isFailure()) return sc;
  sc = m_calCalibSvc->getPed(rngIdxN.getCalXtalId(), pedN, sig, cos);
  if (sc.isFailure()) return sc;

  // convert adc->dac units
  double dacP, dacN;
  sc = m_calCalibSvc->evalDAC(rngIdxP.getCalXtalId(), adcP-pedP, dacP);
  if (sc.isFailure()) return sc;
  sc = m_calCalibSvc->evalDAC(rngIdxN.getCalXtalId(), adcN-pedN, dacN);
  if (sc.isFailure()) return sc;

  // if diodes are different, then I need to convert them both to sm diode
  DiodeNum diodeP = RngNum(rngP).getDiode();
  DiodeNum diodeN = RngNum(rngN).getDiode();
  if (diodeP != diodeN) {
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

      // NOTE:
      // using center of xtal for position since the lrg/sm ratio is 
      // effectively constant throughout the xtal length and this keeps the
      // energy calculation independent of the position calculation.

      double asymNSPB;
      sc = m_calCalibSvc->evalAsymNSPB(xtalId, 0.0, asymNSPB);
      if (sc.isFailure()) return sc;

      dacP /= exp(asymSm - asymNSPB);
      diodeP = SM_DIODE;
    } else { // diodeN == LRG
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

  // retrieve MeVPerDac ratios
  CalibData::ValSig mpdLrg, mpdSm;
  sc = m_calCalibSvc->getMeVPerDac(xtalId, mpdLrg, mpdSm);
  if (sc.isFailure()) return sc;

  double meanDAC = sqrt(dacP*dacN);
  if (diodeP == SM_DIODE)
    energy = meanDAC*mpdSm.getVal();
  else
    energy = meanDAC*mpdLrg.getVal();

  //-- NOISE REDUCTION LAC TEST --/
  // Check that LEX8 ADC values are above 0.5 *'s respective LAC threshold
  belowThresh = false;
  // POS_FACE
  if (rngP == LEX8) {
    CalibData::ValSig fle,fhe,lac;
    // generate tmp xtalid's face specified
    CalXtalId tmpId(xtalId.getTower(),
                    xtalId.getLayer(),
                    xtalId.getColumn(),
                    POS_FACE);
     
    // retreive threshold
    sc = m_calCalibSvc->getTholdCI(tmpId,fle,fhe,lac);

    // set flag
    if (adcP < lac.getVal()*0.5) belowThresh = true;
  }
  // POS_FACE
  if (rngN == LEX8 && !belowThresh) {
    CalibData::ValSig fle,fhe,lac;
    // generate tmp xtalid's face specified
    CalXtalId tmpId(xtalId.getTower(),
                    xtalId.getLayer(),
                    xtalId.getColumn(),
                    NEG_FACE);
     
    // retreive threshold
    sc = m_calCalibSvc->getTholdCI(tmpId,fle,fhe,lac);

    // set flag
    if (adcN < lac.getVal()*0.5) belowThresh = true;
  }

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

  // check that CalXtalId has face & range fields enabled
  if (!xtalId.validFace() || !xtalId.validRange())
    throw invalid_argument("XtalEneTool: CalXtalId requires valid"
                           " face & range fields");

  // retrieve pedestal
  float ped, sig, cos;
  sc = m_calCalibSvc->getPed(xtalId, ped, sig, cos);
  if (sc.isFailure()) return sc;

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

  //-- NOISE REDUCTION LAC TEST --/
  // Check that LEX8 ADC values are above 0.5 *'s respective LAC threshold
  belowThresh = false;
  if (rng == LEX8) {
    CalibData::ValSig fle,fhe,lac;

    // retreive threshold
    sc = m_calCalibSvc->getTholdCI(xtalId,fle,fhe,lac);

    // set flag
    if (adc < lac.getVal()*0.5) belowThresh = true;
  }
  
  return StatusCode::SUCCESS;
}
