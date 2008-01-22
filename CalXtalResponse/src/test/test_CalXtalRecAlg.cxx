// $Header$
/** @file
    @author Z.Fewtrell
*/

// LOCAL INCLUDES
#include "test_CalXtalRecAlg.h"
#include "CalXtalResponse/ICalCalibSvc.h"
#include "CalXtalResponse/ICalSignalTool.h"
#include "../CalCalib/IPrecalcCalibTool.h"
#include "test_util.h"
#include "TestCfg.h"

// GLAST INCLUDES
#include "CalUtil/CalDefs.h"
#include "CalUtil/CalArray.h"
#include "idents/VolumeIdentifier.h"
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"

// EXTLIB INCLUDES
#include "CLHEP/Geometry/Transform3D.h"

// STD INCLUDES
#include <map>
#include <cmath>
#include <algorithm>

using namespace CalUtil;
using namespace std;
using namespace CalXtalResponse;

namespace {
  /// count xtals w/ > 1/2 LAC + margin on both sides 
  short countXtalsOverNoiseThold(IPrecalcCalibTool &precalcCalibTool,
                                 const CalXtalResponse::TwrSet &twrSet,
                                 ICalSignalTool &calSignalTool,
                                 const float marginCIDAC=.2 /**< roughly 2 ADC units @ LEX8**/ ) {
    unsigned short xtalCount = 0;

    /// loop through all active tower bays
    for (TwrSet::const_iterator twrIt(twrSet.begin());
         twrIt != twrSet.end();
         twrIt++) 
      for (tXtalIdx twrXtalIdx; twrXtalIdx.isValid();  twrXtalIdx++) {
        const XtalIdx xtalIdx(*twrIt,twrXtalIdx);

        // check large diode on both faces
        unsigned short nFacesOverThold = 0;
        for (FaceNum face; face.isValid(); face++) {
          const DiodeIdx diodeIdx(xtalIdx, face, LRG_DIODE);
          float diodeSignal;
          if (calSignalTool.getDiodeSignal(diodeIdx, diodeSignal).isFailure())
            return -1;

          const FaceIdx faceIdx(xtalIdx, face);
          float lacCIDAC;
          if (precalcCalibTool.getLacCIDAC(faceIdx, lacCIDAC).isFailure())
            return -1;

          /// quit early if  not over noise threshold on this face
          if (diodeSignal < lacCIDAC/2+marginCIDAC)
            break;

          /// else increment counter
          nFacesOverThold++;
        } // per face loop
        if (nFacesOverThold == FaceNum::N_VALS)
          xtalCount++;

      } // per xtal loop
    return xtalCount;
  }

  //// map XtalIdx to const CalDigi *
  typedef std::map<XtalIdx, const Event::CalDigi *> XtalDigiMap;

  /// return map of xtal Id's to CalDigi *'s
  XtalDigiMap buildXtalDigiMap(const Event::CalDigiCol &calDigiCol) {
    XtalDigiMap retVal;

    for (Event::CalDigiCol::const_iterator digiIter = calDigiCol.begin(); 
         digiIter != calDigiCol.end(); digiIter++)
      retVal[XtalIdx(idents::CalXtalId((*digiIter)->getPackedId()))] = *digiIter;

    return retVal;
  }

  /// check for correct number of readouts in calxtalrecdata (BESTRANGE or ALLRANGE?)
  StatusCode checkNReadouts(const CalXtalResponse::TestCfg &testCfg,
                            const Event::CalXtalRecData &xtalRecData) {
    const unsigned short rightAnswer = (testCfg.trigMode == idents::CalXtalId::BESTRANGE) ? 1 : 4;

    if (xtalRecData.getNReadouts() != rightAnswer)
      return StatusCode::FAILURE;

    return StatusCode::SUCCESS;
  }


  /// calcluate maximum ADC value for given ADC range
  /// \param maxADC destination for max ADC value
  StatusCode maxADC(const CalUtil::RngIdx rngIdx, ICalCalibSvc &calCalibSvc, float &maxADC) {
    // retrieve pedestal calibration
    CalibData::Ped const*const pedCalib = calCalibSvc.getPed(rngIdx);
    if (!pedCalib) return StatusCode::FAILURE;
    const float ped = pedCalib->getAvr();

    vector<float> const*const inlADC = calCalibSvc.getInlAdc(rngIdx);
    if (!inlADC)
      return StatusCode::FAILURE;

    const float maxInlADC = *(inlADC->end()-1);

    maxADC = min<float>(maxInlADC + ped, 4095);

    const RngNum rng = rngIdx.getRng();
    if (rng == HEX1) {
      /// HEX1 saturation marked by ULD value
      const FaceIdx faceIdx(rngIdx.getFaceIdx());
      CalibData::CalTholdCI const*const tholdCICalib = calCalibSvc.getTholdCI(faceIdx);
      if (!tholdCICalib)
        return StatusCode::FAILURE;

      const float uld = tholdCICalib->getULD(rng.val())->getVal();
    
      maxADC = min<float>(uld +ped, maxADC);
    }
    
    return StatusCode::SUCCESS;
  }

  /// possible span of asymmetry across entire crytsal (pos_face_asym - neg_face_asym) (.66 is normal, but can be has high as .9, as low as .6)
  static const float MIN_ASYM_SPAN = .6;
};

/** Algorithm:
    - check correct number of CalXtalRecData (<= # of digis - there
    are possible cuts)
    - check # of recData >= guaranteed numPassCuts
    - verify individual crystal recon
*/
StatusCode test_CalXtalRecAlg::verify(ICalSignalTool &calSignalTool,
                                      ICalCalibSvc &calCalibSvc,
                                      IPrecalcCalibTool &precalcCalibTool,
                                      const CalXtalResponse::TestCfg &testCfg,
                                      const CalXtalResponse::TwrSet &twrSet,
                                      const Event::CalDigiCol &calDigiCol,
                                      const Event::CalXtalRecCol &calXtalRecCol) {

  /// check # of recons <= #of digis
  if (calXtalRecCol.size() > calDigiCol.size()) {
    MsgStream msglog(m_msgSvc, "test_CalXtalRecAlg");   
    msglog << MSG::ERROR << "too many calXtalRecData: " << endreq;
    return StatusCode::FAILURE;
  }

  /// check # of recons >= # of digis certainly over noies suppression thold (lac/2)
  const short xtalsOverNoise = countXtalsOverNoiseThold(precalcCalibTool,
                                                        twrSet,
                                                        calSignalTool);
  if (xtalsOverNoise < 0)
    return StatusCode::FAILURE;

  const XtalDigiMap xtalDigiMap(buildXtalDigiMap(calDigiCol));

  /// verify individual xtals
  for (Event::CalXtalRecCol::const_iterator recIter = calXtalRecCol.begin();
       recIter != calXtalRecCol.end();
       recIter++) {

    const XtalIdx xtalIdx((*recIter)->getPackedId());

    Event::CalDigi const*const calDigi = xtalDigiMap.find(xtalIdx)->second;

    if (verifyXtal(calCalibSvc,
                   testCfg,
                   *calDigi,
                   **recIter).isFailure())
      return StatusCode::FAILURE;
  }

  return StatusCode::SUCCESS;
}


/** Algorithm:
    For each hit xtal:
    - check n readouts.
    - for best range, compare reconstructed values against intended
    deposits in testCfg
    - for all 4 readout ranges (if included) compare recon values
    against digis.
*/
StatusCode test_CalXtalRecAlg::verifyXtal(ICalCalibSvc &calCalibSvc,
                                          const CalXtalResponse::TestCfg &testCfg,
                                          const Event::CalDigi &calDigi,
                                          const Event::CalXtalRecData &calXtalRecData) {
  
  if (checkNReadouts(testCfg, calXtalRecData).isFailure())
    return StatusCode::FAILURE;

  const XtalIdx xtalIdx(idents::CalXtalId(calXtalRecData.getPackedId()));

  /// check best range against MC truth (other ranges just check against associated 
  /// digi
  if (verifyRangeRecon(xtalIdx,
                       calCalibSvc,
                       testCfg,
                       *calDigi.getXtalReadout(0),
                       *calXtalRecData.getRangeRecData(0)).isFailure())
    return StatusCode::FAILURE;

  // check each individual readout (check that recon matches digi for each redout)
  for (unsigned short nReadout = 0;
       nReadout < calXtalRecData.getNReadouts();
       nReadout++) {
    if (verifyRangeReadout(xtalIdx,
                           calCalibSvc,
                           *calDigi.getXtalReadout(nReadout),
                           *calXtalRecData.getRangeRecData(nReadout)).isFailure())
      return StatusCode::FAILURE;
  }
  
  return StatusCode::SUCCESS;
}

/** Algorithm:
    From Digi data:
    - retrieve ADC & range info from digi
    - check that ADC value is <= saturation value
    - pedestal subtract adc values
    - quit early with SUCCESS if adc < 1 (we use sqrt & log() &
    integer value of adc=0 precludes this.
    - evaluate faceSignal (MeV) for each xtal face
    - evaluate light asymmetry from faceSignal measurements
    - evaluate deposit energy from geom. mean of faceSignal
    measurements
    From Recon Data:
    - retrieve meV & position & compare against values determined from
    Digi.

    - allow large tolerance for hits involving HE diode signal as
    ground-muon based asymmetry calculations on HE diode have high
    error bars.
    - for position comparison, convert percent tolerance to mm along
    crystal end.
*/
StatusCode test_CalXtalRecAlg::verifyRangeReadout(const CalUtil::XtalIdx xtalIdx,
                                                  ICalCalibSvc &calCalibSvc,
                                                  const Event::CalDigi::CalXtalReadout &digiRO,
                                                  const Event::CalXtalRecData::CalRangeRecData &reconRO) {

  /// which diode readout was used on each face?
  CalArray<FaceNum, DiodeNum> diode;
  /// what was the CIDAC signal level on that face
  CalArray<FaceNum, float> cidac;
  CalArray<FaceNum, float> faceSignal;
  /// pedestal subtracted adc for each face
  CalArray<FaceNum, float> adcPed;

  for (FaceNum face; face.isValid();  face++) {
    // get readout range number & adc values
    const RngNum rng = reconRO.getRange(face);
    const RngIdx rngIdx(xtalIdx, face, rng);
    diode[face] = rng.getDiode();
    const float adc = digiRO.getAdc(face);

    /// check for saturated ADC channels.
    float max_adc;
    if (maxADC(rngIdx, calCalibSvc, max_adc).isFailure())
      return StatusCode::FAILURE;
    if (adc >= floor(max_adc))
      return StatusCode::SUCCESS;

    // retrieve pedestal calibration
    CalibData::Ped const*const pedCalib = calCalibSvc.getPed(rngIdx);
    if (!pedCalib) return 0;
    const float ped = pedCalib->getAvr();

    // ped subtracted ADC
    adcPed[face] = adc - ped;

    // skip signals < 1 adc unit
    if (adcPed[face] < 1)
      return StatusCode::SUCCESS;

    if (calCalibSvc.evalCIDAC(rngIdx, adcPed[face], cidac[face]).isFailure())
      return StatusCode::FAILURE;

    if (calCalibSvc.evalFaceSignal(rngIdx, adcPed[face], faceSignal[face]).isFailure())
      return StatusCode::FAILURE;    
  }

  const float testAsym = log(cidac[POS_FACE]/cidac[NEG_FACE]);
  /// reconstructed position in mm from crystal center.
  float testPosMM;
  if (calCalibSvc.evalPos(xtalIdx, 
                          AsymType(diode[POS_FACE],diode[NEG_FACE]),
                          testAsym,
                          testPosMM).isFailure())
    return StatusCode::FAILURE;

  // 1st clip position to end of xtal
  testPosMM = min<float>(   m_csiLength/2,testPosMM);
  testPosMM = max<float>(-1*m_csiLength/2,testPosMM);


  const float testEne = sqrt(faceSignal[POS_FACE]*faceSignal[NEG_FACE]);
  
  /// test energy (value should be within 1 adc unit of real answer)
  const float recEne = reconRO.getEnergy();

  /// mixed diode results can suffer from incosistent HE asymmetry curves (due to muon calibration)
  const float pctTolerance = (diode[POS_FACE] == diode[NEG_FACE]) ?
    MAX_RECON_DIFF : MAX_ASYM_DIFF;

  if (!smart_compare(recEne, testEne, pctTolerance)) {
    MsgStream msglog(m_msgSvc, "test_CalXtalRecAlg");   
    msglog << MSG::ERROR << "bad reconstructed energy: " << xtalIdx.toStr() << endreq;
    return StatusCode::FAILURE;
  }

  const Point recPosLAT = reconRO.getPosition();
  const Point testPosLAT = pos2Point(xtalIdx, testPosMM);
  const float posDiffMM = (recPosLAT - testPosLAT).magnitude();

  /// convert tolerance from pct of asymmetry to units of one crystal length.
  static const float posToleranceXtal = pctTolerance/MIN_ASYM_SPAN;
  const float posToleranceMM = posToleranceXtal*m_csiLength;

  if (posDiffMM > posToleranceMM) {
    MsgStream msglog(m_msgSvc, "test_CalXtalRecAlg");   
    msglog << MSG::ERROR << "bad reconstructed position: " << xtalIdx.toStr() << endreq;
    return StatusCode::FAILURE;
  }

  return StatusCode::SUCCESS;
}


/** Algorithm:
    - check for digi saturation & quit early with success is we are
    saturated.
    - pedestal subtract adc values
    - retrieve intended deposit energy & position from testCfg
    - retrieve reconstructed deposit energy & position for
    CalRangeRecData
    - deterimine quantization error from floor() of pedestal
    subtracted ADC
    - update tolerance based on quantization error
    - allow large tolerance for hits involving HE diode signal as
    ground-muon based asymmetry calculations on HE diode have high
    error bars.
    - for position comparison, convert percent tolerance to mm along
    crystal end.
    
*/
StatusCode test_CalXtalRecAlg::verifyRangeRecon(const CalUtil::XtalIdx xtalIdx,
                                                ICalCalibSvc &calCalibSvc,
                                                const CalXtalResponse::TestCfg &testCfg,
                                                const Event::CalDigi::CalXtalReadout &digiRO,
                                                const Event::CalXtalRecData::CalRangeRecData &reconRO) {
  
  
  /// keep track of diode used on each xtal face
  CalArray<FaceNum, DiodeNum> diode;
  /// pedestal subtracted adc on each xtal face
  CalArray<FaceNum, float> adcPed;

  /// check for digi saturation
  for (FaceNum face; face.isValid();  face++) {
    const unsigned short adc = digiRO.getAdc(face);
    const RngNum rng(digiRO.getRange(face));
    const RngIdx rngIdx(xtalIdx, face,rng);
    diode[face] = rng.getDiode();

    float max_adc;
    if (maxADC(rngIdx, calCalibSvc, max_adc).isFailure())
      return StatusCode::FAILURE;

    if (adc >= floor(max_adc))
      return StatusCode::SUCCESS;

    // retrieve pedestal calibration
    CalibData::Ped const*const pedCalib = calCalibSvc.getPed(rngIdx);
    if (!pedCalib) return StatusCode::FAILURE;
    const float ped = pedCalib->getAvr();

    adcPed[face] = adc - ped;
  }

  /// compare to simulation intent.
  const float intentEne = testCfg.totalEnergy();

  /// intended position along length of xtal
  const float intentWeightedPosRel = testCfg.weightedPos();
  /// in mm from center of xtal
  const float intentWeightedPosMM = (intentWeightedPosRel - .5)*m_csiLength;
  /// intended posiiton in LAT coordinates.
  const Point intentPos = pos2Point(xtalIdx, intentWeightedPosMM);

  /// retrieve recon info
  const float recEne = reconRO.getEnergy();
  const Point recPos = reconRO.getPosition();

  // quantization error
  const float minAdcPed = floor(min(adcPed[POS_FACE], adcPed[NEG_FACE]));
  const float quantError = 1/minAdcPed;
  float pctTolerance = max(MAX_RECON_DIFF, quantError);

  // curvature in asymmetry model introduces error in summed energy deposits.
  if (testCfg.xtalHits.size() > 1)
    pctTolerance = max(pctTolerance, MAX_ASYM_DIFF);

  /// mixed diode results can suffer from incosistent HE asymmetry curves (due to muon calibration)
  if (diode[POS_FACE] != diode[NEG_FACE])
    pctTolerance = max(pctTolerance, MAX_ASYM_DIFF);

  /// HE diode asymmetry can be pretty bad when measured with muons.
  if (find(diode.begin(), diode.end(),SM_DIODE) != diode.end())
    pctTolerance = max(pctTolerance, MAX_ASYM_DIFF);


  if (!smart_compare(recEne, intentEne, pctTolerance)) {
    MsgStream msglog(m_msgSvc, "test_CalXtalRecAlg");   
    msglog << MSG::ERROR << "bad reconstructed energy: " << xtalIdx.toStr() << endreq;
    return StatusCode::FAILURE;
  }

  /// get position difference in mm
  const float posDiff = Vector(recPos - intentPos).magnitude();

  /// tolerance in units of one crystal length.
  /// asymmetry ratio on one side is roughly .65 > than the other side.
  const float posToleranceXtal = pctTolerance/MIN_ASYM_SPAN;
  const float posToleranceMM = posToleranceXtal*m_csiLength;
  
  if (posDiff > posToleranceMM) {
    MsgStream msglog(m_msgSvc, "test_CalXtalRecAlg");   
    msglog << MSG::ERROR << "bad reconstructed position: " << xtalIdx.toStr() << endreq;
    return StatusCode::FAILURE;
    
  }

  return StatusCode::SUCCESS;
}


Point test_CalXtalRecAlg::pos2Point(const XtalIdx xtalIdx, 
                                    const float xtalPosMM) const {
  //-- CONSTRUCT GEOMETRY VECTORS FOR XTAL --//
    
  // create Volume Identifier for segments 0 & 11 of this crystal
    
  // volId info snagged from 
  // http://www.slac.stanford.edu/exp/glast/ground/software/geometry/docs/identifiers/geoId-RitzId.shtml
  idents::VolumeIdentifier segm0Id, segm11Id;

  const TwrNum twr = xtalIdx.getTwr();
  const LyrNum lyr = xtalIdx.getLyr();
  const ColNum col = xtalIdx.getCol();
  
  // init seg0 w/ info shared by both.
  segm0Id.append(m_eLATTowers);
  segm0Id.append(twr.getRow());
  segm0Id.append(twr.getCol());
  segm0Id.append(m_eTowerCAL);
  segm0Id.append(lyr.val());
  segm0Id.append(lyr.getDir().val()); 
  segm0Id.append(col.val());
  segm0Id.append(m_eXtal);

  // copy over shared info
  segm11Id = segm0Id;
  
  // init segment specific info.
  segm0Id.append(0); // segment Id
  segm11Id.append(m_nCsISeg-1); // set segment number for the last
  // segment
  
  HepTransform3D transf;

  //get 3D transformation for segment 0 of this crystal
  m_detSvc->getTransform3DByID(segm0Id,&transf);
  //get position of the center of the segment 0
  Vector vect0 = transf.getTranslation();

  //get 3D transformation for the last segment of this crystal
  m_detSvc->getTransform3DByID(segm11Id,&transf);
  //get position of the center of the last segment
  Vector vect11 = transf.getTranslation();

  Point p0(0.,0.,0.);           
  // position of the crystal center
  Point pCenter = p0+(vect0+vect11)*0.5; 
  //normalized vector of the crystal direction 
  Vector dirXtal = (vect11-vect0)*m_nCsISeg/(m_nCsISeg-1);  

  // put 1D position info into 3D vector
  // 'xtalPos' is in units of xtal Length, convert to rel units
  // (-1->1)
  const float xtalPosRel = xtalPosMM / m_csiLength; 
  return pCenter+dirXtal*xtalPosRel;

}

/// Gaudi initialize() phase intializations
StatusCode test_CalXtalRecAlg::initialize(IMessageSvc *msgSvc,
                                          IGlastDetSvc *detSvc) {
  m_msgSvc = msgSvc;
  m_detSvc = detSvc;

  //-- RETRIEVE CONSTANTS --//

  double tmp;
  // map containing pointers to integer constants to be read
  // with their symbolic names from xml file used as a key 
  typedef map<int*,string> PARAMAP;
  PARAMAP param;

  // INT CONSTANTS
  //     filling the map with information on constants to be read 
  param[&m_eTowerCAL]    = string("eTowerCAL");
  param[&m_eLATTowers]   = string("eLATTowers");
  param[&m_nCsISeg]      = string("nCsISeg");
  param[&m_eXtal]        = string("eXtal");

  // loop over all constants information contained in the map
  for(PARAMAP::iterator iter=param.begin(); iter!=param.end();iter++){
    //  retrieve constant
    if(!m_detSvc->getNumericConstByName((*iter).second, &tmp)) {
      MsgStream msglog(m_msgSvc, "test_CalXtalRecAlg");
      msglog << MSG::ERROR << " constant " <<(*iter).second
             << " not defined" << endreq;
      return StatusCode::FAILURE;
    } else *((*iter).first)= int(tmp); // store retrieved value 
  }

  // DOUBLE CONSTANTS
  if (m_detSvc->getNumericConstByName("CsILength", &tmp).isFailure()) {
    MsgStream msglog(m_msgSvc, "test_CalXtalRecAlg");   
    msglog << MSG::ERROR << " constant CsILength not defined" << endreq;
    return StatusCode::FAILURE;
  }
  m_csiLength = tmp;

  return StatusCode::SUCCESS;
}
