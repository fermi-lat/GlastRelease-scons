// $Header$
/** @file
    @author Z.Fewtrell
*/

// LOCAL
#include "XtalRecTool.h"
#include "CalXtalResponse/ICalCalibSvc.h"
#include "../Xtalk/INeighborXtalkTool.h"

// GLAST
#include "idents/VolumeIdentifier.h"
#include "geometry/Point.h"
#include "CalUtil/CalVec.h"
#include "CalUtil/CalArray.h"
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"


#include "Event/Digi/CalDigi.h"
#include "Event/TopLevel/Event.h"


// EXTLIB
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SmartDataPtr.h"

#include "CLHEP/Geometry/Transform3D.h"

// STD
#include <cmath>
#include <map>
#include <algorithm>
#include <memory>

using namespace CalUtil;
using namespace Event;
using namespace CalibData;
using namespace idents;
using namespace std;

//static ToolFactory<XtalRecTool> s_factory;
//const IToolFactory& XtalRecToolFactory = s_factory;
DECLARE_TOOL_FACTORY(XtalRecTool);

XtalRecTool::XtalRecTool( const string& type,
                          const string& name,
                          const IInterface* parent)
  : AlgTool(type,name,parent),
    m_calCalibSvc(0),
    m_CsILength(0),
    m_eLATTowers(-1),
    m_eTowerCAL(-1),
    m_eXtal(-1),
    m_nCsISeg(-1)
{
  declareInterface<IXtalRecTool>(this);

  declareProperty("CalCalibSvc", m_calCalibSvcName="CalCalibSvc");
}

/// init / retrieve all needed Gaudi objects
StatusCode XtalRecTool::initialize() {
  MsgStream msglog(msgSvc(), name());
  msglog << MSG::INFO << "initialize" << endreq;

  StatusCode sc;

  //-- jobOptions --//
  if ((sc = setProperties()).isFailure()) {
    msglog << MSG::ERROR << "Failed to set properties" << endreq;
    return sc;
  }
  
  // obtain CalCalibSvc
  sc = service(m_calCalibSvcName.value(), m_calCalibSvc);
  if (sc.isFailure()) {
    msglog << MSG::ERROR << "can't get " << m_calCalibSvcName  << endreq;
    return sc;
  }

  //-- EXTRACT INT CONSTANTS --//
  double tmp;
  // map containing pointers to integer constants to be read
  // with their symbolic names from xml file used as a key 
  typedef map<int*,string> PARAMAP;
  PARAMAP param;

  //     filling the map with information on constants to be read 
  param[&m_eTowerCAL]    = string("eTowerCAL");
  param[&m_eLATTowers]   = string("eLATTowers");
  param[&m_nCsISeg]      = string("nCsISeg");
  param[&m_eXtal]        = string("eXtal");
    
  //-- RETRIEVE GlastDevSvc --//

  sc = service("GlastDetSvc", m_detSvc);
  // loop over all constants information contained in the map
  for(PARAMAP::iterator iter=param.begin(); iter!=param.end();iter++){
    //  retrieve constant
    if(!m_detSvc->getNumericConstByName((*iter).second, &tmp)) {
      MsgStream msglog(msgSvc(), name());
      msglog << MSG::ERROR << " constant " <<(*iter).second
             <<" not defined" << endreq;
      return StatusCode::FAILURE;
    } else *((*iter).first)= int(tmp); // store retrieved value 
  }

  //-- EXTRACT DOUBLE CONSTANTS --//

  // doubles are done separately
  sc = m_detSvc->getNumericConstByName("CsILength", &tmp);
  if (sc.isFailure()) {
    msglog << MSG::ERROR << " constant CsILength not defined" << endreq;
    return sc;
  }
  m_CsILength = tmp;

  return StatusCode::SUCCESS;
}

/**
   - convert adc readouts to CIDAC scale
   - compute centroid position & intesity from CIDAC scale
   - check for noise threshold & xtal saturation
*/
StatusCode XtalRecTool::calculate(const Event::CalDigi &         digi,
                                  Event::CalXtalRecData&         xtalRec,
                                  CalVec<FaceNum, bool>&         belowNoise,
                                  CalVec<FaceNum, bool>&         saturated,
                                  CalVec<FaceNum, bool>&         adcSaturated,
                                  INeighborXtalkTool const*const xtalkTool) 
{
    // initialize return values
    fill(belowNoise.begin(), belowNoise.end(), false);
    fill(saturated.begin(), saturated.end(), false);
    fill(adcSaturated.begin(), adcSaturated.end(), false);
  
    // used for global access routines.
    const XtalIdx xtalIdx(digi.getPackedId());

    // check for empty digi
    if (digi.getReadoutCol().size() <= 0)
        return StatusCode::SUCCESS;

    for (CalDigi::CalXtalReadoutCol::const_iterator ro =  digi.getReadoutCol().begin();
         ro != digi.getReadoutCol().end();
         ro++) 
    {
        /// check if current range is below noise threshold
        CalVec<FaceNum, bool> rngBelowNoise;

        // using auto_ptr, bc I am responsible for this pointer.
        auto_ptr<Event::CalXtalRecData::CalRangeRecData> rngRec(createRangeRecon(xtalIdx,
                                                                                 *ro, 
                                                                                 rngBelowNoise,
                                                                                 saturated,
                                                                                 adcSaturated,
                                                                                 xtalkTool));
    
        if (!rngRec.get())
            return StatusCode::FAILURE;

        /// only retain belowNoise info for first range
        if (ro == digi.getReadoutCol().begin()){
          copy(rngBelowNoise.begin(), rngBelowNoise.end(), belowNoise.begin());
        }

        xtalRec.addRangeRecData(*(rngRec.get()));

        /// Finally: propagate the ADC saturation information.
        /// We've been discussing whether we should use the saturated variable, here, and the
        /// resolution was that we shouldn't. Decision was taken to pick a value slightly lower
        /// than the 12-bit ADC maximum value, i.e. 4060 to be consistent with the value
        /// historically choosen by Philippe Bruel in the original implementation.
        /// That's the reason for introducing the new variable adcSaturated.
        /// Luca Baldini, November 27, 2012.
        if (adcSaturated[POS_FACE]){
          xtalRec.setPosFaceSaturated();
        }
        if (adcSaturated[NEG_FACE]){
          xtalRec.setNegFaceSaturated();
        }
        
    } // per range loop

    return StatusCode::SUCCESS;
}

// Algorithm to deal with Ambiguous Region of Xtal response
StatusCode XtalRecTool::ambiguity(Event::CalXtalRecCol* pxtalrecs)
{
    if (!pxtalrecs->empty()) 
    {
        for(Event::CalXtalRecCol::const_iterator jlog=pxtalrecs->begin(); jlog != pxtalrecs->end(); ++jlog) 
        {
            Event::CalXtalRecData& recLog = **jlog;   
            idents::CalXtalId xtalId = recLog.getPackedId();
            int tower = xtalId.getTower();
            int layer = xtalId.getLayer();
            Point pos = recLog.getPosition();
            // Hardwire in is that even layers have the long. co-ordinate in the x direction....
            float xLong;
			point2Pos(xtalId, xLong, pos);
            float xTrans; 
            
            double eneMax = 0.; 

            for(Event::CalXtalRecCol::const_iterator klog=pxtalrecs->begin(); klog != pxtalrecs->end(); ++klog) 
            {
                Event::CalXtalRecData& recLogTs = **klog;   
                idents::CalXtalId xtalIdTs = recLogTs.getPackedId();
                int towerTs = xtalIdTs.getTower();
                if(towerTs != tower) continue;

                int layerTs = xtalIdTs.getLayer();
                if(layerTs == layer+1 || layerTs == layer-1) 
                {
                    double eneLog = recLogTs.getEnergy();
                    if(eneLog > eneMax) 
                    {
                        eneMax = eneLog;
                        Point posTs = recLogTs.getPosition();
						point2PosTrans(xtalIdTs, xTrans, posTs);
                    }
                }
            }
			if(eneMax > 5. && fabs(xTrans) > 135. && fabs(xLong) < 135.) {
                double deltaX = 135. - fabs(xLong);
                double corrTerm = deltaX*deltaX*.15;
                float xLongCorr = (xTrans > 0) ? xLong + corrTerm : xLong - corrTerm;

				// Clip to max end of xtal
				xLongCorr = min<float>(    m_CsILength/2, xLongCorr);
                xLongCorr = max<float>(-1.*m_CsILength/2, xLongCorr);

                Point corrPos(0., 0., 0.);
				const CalUtil::XtalIdx xtalIdxU(xtalId);
				pos2Point(xtalIdxU, xLongCorr, corrPos);

                const Event::CalXtalRecData::CalRangeRecData *rngData = recLog.getRangeRecData(0); 
				Event::CalXtalRecData::CalRangeRecData *rngDataR = const_cast<Event::CalXtalRecData::CalRangeRecData *> (rngData);
                rngDataR->setPosition(corrPos);
            }
        }
    }
    return StatusCode::SUCCESS;
}


/*   Algorithm is as follows:
   - subtract pedestals from ADC readouts
   - check if either face is below 1/2 LAC signal level
   - convert ADC readouts to CIDAC scale
   - detect saturated HEX1 range
   - (optional) neighboring crytal cross-talk correction
   - use asymmetry ratio of both cidac signal to evaluate longitudinal
   deposit centroid.
   - use geometric mean of both cidac signals to evaluate deposit
   energy.
   - create CalRangeRecData TDS object
 */
Event::CalXtalRecData::CalRangeRecData *XtalRecTool::createRangeRecon(const CalUtil::XtalIdx xtalIdx,
                                                                      const Event::CalDigi::CalXtalReadout &ro,
                                                                      CalVec<FaceNum, bool> &belowNoise,
                                                                      CalVec<FaceNum, bool> &saturated,
                                                                      CalVec<FaceNum, bool> &adcSaturated,
                                                                      INeighborXtalkTool const*const xtalkTool
                                                                      ) const {
  StatusCode sc;

  CalArray<FaceNum, DiodeNum> diode;
  CalArray<FaceNum, float> cidac;
  CalArray<FaceNum, RngNum> rng;
  CalArray<FaceNum, float> adcPed;
  
  //////////////////////////////////////
  //-- STEP 1: PEDESTAL SUBTRACTION --//
  //////////////////////////////////////
  for (FaceNum face; face.isValid();  face++) {
    // get readout range number & adc values
    rng[face] = ro.getRange(face);
    const float adc = ro.getAdc(face);
    const RngIdx rngIdx(xtalIdx, face, rng[face]);
    diode[face] = rng[face].getDiode();

    // retrieve pedestal calibration
    float ped,pedSig;
    StatusCode sc = m_calCalibSvc->getPed(rngIdx,ped);
    if (sc.isFailure()) return 0;
    sc = m_calCalibSvc->getPedSig(rngIdx,pedSig);
    if (sc.isFailure()) return 0;

    // ped subtracted ADC
    adcPed[face] = adc - ped;
  
    /////////////////////////////////
    //-- STEP 2: NOISE REDUCTION --//
    /////////////////////////////////
    const FaceIdx faceIdx(xtalIdx, face);
    CalTholdCI const*const tholdCI = m_calCalibSvc->getTholdCI(faceIdx);
    // get LAC threshold calib
    if (!tholdCI) return 0;
    if (rng[face] == LEX8) {
      const float lacThresh = tholdCI->getLAC()->getVal();
    
      // LEX8 range is compared against 0.5 * lac threshold
      // we throw out entire xtal if adc is too low 
      if (adcPed[face] < lacThresh*0.5)
        belowNoise[face] = true;
    } else {
      // other ranges are compared against 5 sigma energy 
      if (adcPed[face] < pedSig*5.0)
        belowNoise[face] = true;
    }

    /////////////////////////////////////////////
    //-- STEP 3: CONVERT ADCs -> CIDAC SCALE --//
    /////////////////////////////////////////////
    sc = m_calCalibSvc->evalCIDAC(rngIdx, adcPed[face], cidac[face]);
    if (sc.isFailure()) return 0;
        
    // check for invalid cidac vals (i need to take the sqrt AND the log)
    if (cidac[face] <= 0) {
      // create MsgStream only when needed (for performance)
      MsgStream msglog(msgSvc(), name()); 
      msglog << MSG::VERBOSE;
      // need to use .stream() to get xtalId to pretty print.
      msglog.stream() << "minimal signal... CIDAC val <= 0, can't calculate energy/position." 
                      <<  " xtal=[" << xtalIdx.getCalXtalId() << ']';
      msglog << endreq;

      belowNoise[face] = true;

      // this is not a FAILURE condition, just not enough signal
      // skip to next range
      continue;
    }

    /////////////////////////////////////////////
    //-- STEP 4: DETECT SATURATED HEX1 RANGE --//
    /////////////////////////////////////////////
    if (rng[face] == HEX1) {
      const float h1Limit = tholdCI->getULD(HEX1.val())->getVal();
      if (adc >= h1Limit) saturated[face] = true;
      if (adc >= CAL_ADC_SATURATION) adcSaturated[face] = true;
    }

    ///////////////////////////////////////////////////////////
    //-- STEP 5: (OPTIONAL) NEIGHBOR XTAL XTALK CORRECTION --//
    ///////////////////////////////////////////////////////////
    if (xtalkTool 
        && !saturated[face])  // don't subtract signal from saturated xtals
      {
        const DiodeIdx diodeIdx(xtalIdx, face, diode[face]);
        float xtalkCIDAC;
        sc = xtalkTool->calcXtalkCIDAC(diodeIdx, xtalkCIDAC);
        if (sc.isFailure()) return 0;
        cidac[face] -= xtalkCIDAC;
      }
  }

  ////////////////////////////////////
  //-- STEP 7: CALCULATE POSITION --//
  ////////////////////////////////////

  const float asym = log(cidac[POS_FACE]/cidac[NEG_FACE]);
  float pos;
  if (diode[POS_FACE] == LRG_DIODE)
    if (diode[NEG_FACE] == LRG_DIODE)
      sc = m_calCalibSvc->evalPos(xtalIdx, ASYM_LL, asym, pos);   
    else
      sc = m_calCalibSvc->evalPos(xtalIdx, ASYM_LS, asym, pos);
  else
    if (diode[NEG_FACE] == SM_DIODE)
      sc = m_calCalibSvc->evalPos(xtalIdx, ASYM_SS, asym, pos);
    else
      sc = m_calCalibSvc->evalPos(xtalIdx, ASYM_SL, asym, pos);
  if (sc.isFailure()) return 0;

  // First attempts to correct overshoot 
  if(fabs(pos) > 113.) {
    float deltaX = fabs(pos)-113.;
    float deltaCorrection = deltaX*deltaX*.0052;
	if(fabs(pos) - deltaCorrection < 113) deltaX = 0.;  //Stop unwanted large deltaX behavior
    if(pos > 0. ) pos -= deltaCorrection;
    else          pos += deltaCorrection;
  }
  // Now clip position to end of xtal
    pos = min<float>(   m_CsILength/2,pos);
    pos = max<float>(-1.*m_CsILength/2,pos);

  
  // generate position vector from scalar longitudinal position
  Point pXtal;
  pos2Point(xtalIdx, pos, pXtal);


  ///////////////////////////////////////
  //-- STEP 8: CALCULATE TRUE ENERGY --//
  ///////////////////////////////////////
  
  //-- convert diodes to same size (if needed)
  DiodeNum mpdDiode(diode[POS_FACE]); // diode size to use for MeVPerDAC conversion
  CalVec<FaceNum,float> mpdCIDAC;
  mpdCIDAC[POS_FACE]  = cidac[POS_FACE];   // cidac value to use for MeVPerCIDAC conversion
  mpdCIDAC[NEG_FACE]  = cidac[NEG_FACE];   // cidac value to use for MeVPerCIDAC conversion

  // if diodes on each face differ convert to small on both faces
  if (diode[POS_FACE].val() != diode[NEG_FACE].val()) {
    if (diode[POS_FACE] == LRG_DIODE)
      sc = largeCIDAC2Small(xtalIdx,
                            POS_FACE, 
                            pos, 
                            mpdCIDAC[POS_FACE], 
                            mpdCIDAC[POS_FACE]);
    else // diode[NEG_FACE] == LRG_DIODE
      sc = largeCIDAC2Small(xtalIdx,
                            NEG_FACE, 
                            pos, 
                            mpdCIDAC[NEG_FACE], 
                            mpdCIDAC[NEG_FACE]);
    
    if (sc.isFailure()) return 0;
    
    mpdDiode = SM_DIODE;
  }

  //-- RETRIEVE MEV PER DAC--// 
  CalMevPerDac const*const mpdCalib = m_calCalibSvc->getMPD(xtalIdx);
  if (!mpdCalib) return 0;
  CalVec<DiodeNum, float> mpd;
  mpd[LRG_DIODE] = mpdCalib->getBig()->getVal();
  mpd[SM_DIODE]  = mpdCalib->getSmall()->getVal();


  const float meanCIDAC = sqrt(mpdCIDAC[POS_FACE]*mpdCIDAC[NEG_FACE]);
  const float ene = meanCIDAC*mpd[mpdDiode];

  // just to try
  float eneP = mpd[mpdDiode] * mpdCIDAC[POS_FACE];
  float eneN = mpd[mpdDiode] * mpdCIDAC[NEG_FACE];

  ////////////////////////////////////
  //-- STEP 10: POPULATE TDS CLASS --//
  ////////////////////////////////////
  auto_ptr<CalXtalRecData::CalRangeRecData> rngRec(new CalXtalRecData::CalRangeRecData((CalXtalId::AdcRange)rng[POS_FACE], eneP, 
                                                                                       (CalXtalId::AdcRange)rng[NEG_FACE], eneN));
  rngRec->setPosition(pXtal);
  
  return rngRec.release();
}

StatusCode XtalRecTool::largeCIDAC2Small(const XtalIdx xtalIdx,
                                         const FaceNum face, 
                                         const float pos, 
                                         const float largeCIDAC, 
                                         float &smallCIDAC) const {
  StatusCode sc;

  // small diode asymmetry is used for both if() cases
  float asymSm;
  sc = m_calCalibSvc->evalAsym(xtalIdx, ASYM_SS, pos, asymSm);
  if (sc.isFailure()) return sc;


  if (face == POS_FACE) {
    // STRATEGY:
    // convert PL -> PS
    // AsymNSPB = log(PL / NS)
    // AsymSm   = log(PS / NS)
    // exp(asymSm)/exp(asymNSPB) = (PS/NS)/(PL/NS) = PS/PL
    //      exp(asymSm-asymNSPB) = PS/PL
    //                        PS = PL*exp(asymSm=asymNSPB
      
    float asymNSPB;    
    sc = m_calCalibSvc->evalAsym(xtalIdx, ASYM_LS, pos, asymNSPB);
    if (sc.isFailure()) return sc;

    smallCIDAC = largeCIDAC * exp(asymSm - asymNSPB);
  } 
  else { 
    // STRATEGY:
    // convert NL -> NS
    // asymPSNB = log(PS/NL)
    // asymSm   = log(PS/NS)
    // exp(asymPSNB)/exp(asymSm) = (PS/NL)/(PS/NS)
    //      exp(asymPSNB-asymSm) = NS/NL
    //                        NS = NL*exp(asymPSNB-asymSm)

    float asymPSNB;
    sc = m_calCalibSvc->evalAsym(xtalIdx, ASYM_SL, pos, asymPSNB);
    if (sc.isFailure()) return sc;

    smallCIDAC = largeCIDAC * exp(asymPSNB - asymSm);
  }

  return StatusCode::SUCCESS;

}


/** \brief convert longitudinal pos along cal xtal to 3d point in LAT geometry space

    \param pXtal ouput position vector
    \param pos input longitudinal position in mm from center of xtal
*/
void XtalRecTool::pos2Point(const XtalIdx xtalIdx,
                            const float pos, 
                            Point &pXtal) const {
  //-- CONSTRUCT GEOMETRY VECTORS FOR XTAL --//

  // create Volume Identifier for segments 0 & 11 of this crystal
  // volId info snagged from 
  // http://www.slac.stanford.edu/exp/glast/ground/software/geometry/docs/identifiers/geoId-RitzId.shtml
  idents::VolumeIdentifier segm0Id, segm11Id;
  
  // init seg0 w/ info shared by both.
  segm0Id.append(m_eLATTowers);
  segm0Id.append(xtalIdx.getTwr().getRow());
  segm0Id.append(xtalIdx.getTwr().getCol());
  segm0Id.append(m_eTowerCAL);
  segm0Id.append(xtalIdx.getLyr().val());
  segm0Id.append(xtalIdx.getLyr().getDir().val()); 
  segm0Id.append(xtalIdx.getCol().val());
  segm0Id.append(m_eXtal);

  // copy over shared info
  segm11Id = segm0Id;
  
  // init segment specific info.
  segm0Id.append(0); // segment Id
  segm11Id.append(m_nCsISeg-1); // set segment number for the last segment
  
  HepTransform3D transf;

  //get 3D transformation for segment 0 of this crystal
  m_detSvc->getTransform3DByID(segm0Id,&transf);
  //get position of the center of the segment 0
  const Vector vect0 = transf.getTranslation();

  //get 3D transformation for the last segment of this crystal
  m_detSvc->getTransform3DByID(segm11Id,&transf);
  //get position of the center of the last segment
  const Vector vect11 = transf.getTranslation();
          
  Point p0(0., 0., 0.);
  // position of the crystal center - need to bogus p0 to cast Vector -> Point
  const Point pCenter = p0 + (vect0+vect11)*0.5; 
  //normalized vector of the crystal direction 
  const Vector dirXtal = (vect11-vect0).unit();  

  // put 1D position info into 3D vector
  pXtal = pCenter+dirXtal*pos;
}

/** \brief convert 3d point in LAT geometry space to longitudinal pos along cal xtal
    \param pXtal ouput position vector
    \param pos input longitudinal position in mm from center of xtal
*/
void XtalRecTool::point2Pos(const idents::CalXtalId xtalId,
                            float &pos, 
                            const Point pXtal) const {
  
  // create Volume Identifier for segments 0 & 11 of this crystal
  // volId info snagged from 
  // http://www.slac.stanford.edu/exp/glast/ground/software/geometry/docs/identifiers/geoId-RitzId.shtml

  const CalUtil::XtalIdx xtalIdx(xtalId);
  idents::VolumeIdentifier segm0Id, segm11Id;
  
  // init seg0 w/ info shared by both.
  segm0Id.append(m_eLATTowers);
  segm0Id.append(xtalIdx.getTwr().getRow());
  segm0Id.append(xtalIdx.getTwr().getCol());
  segm0Id.append(m_eTowerCAL);
  segm0Id.append(xtalIdx.getLyr().val());
  segm0Id.append(xtalIdx.getLyr().getDir().val()); 
  segm0Id.append(xtalIdx.getCol().val());
  segm0Id.append(m_eXtal);

  // copy over shared info
  segm11Id = segm0Id;
  
  // init segment specific info.
  segm0Id.append(0); // segment Id
  segm11Id.append(m_nCsISeg-1); // set segment number for the last segment
  
  HepTransform3D transf;

  //get 3D transformation for segment 0 of this crystal
  m_detSvc->getTransform3DByID(segm0Id,&transf);
  //get position of the center of the segment 0
  const Vector vect0 = transf.getTranslation();

  //get 3D transformation for the last segment of this crystal
  m_detSvc->getTransform3DByID(segm11Id,&transf);
  //get position of the center of the last segment
  const Vector vect11 = transf.getTranslation();
          
  Point p0(0., 0., 0.);
  // position of the crystal center - need to bogus p0 to cast Vector -> Point
  const Point pCenter = p0 + (vect0+vect11)*0.5; 

  Point hit = p0 +(pXtal - pCenter);
  pos = (xtalId.isX()) ?  hit.x() : hit.y();
}

void XtalRecTool::point2PosTrans(const idents::CalXtalId xtalId,
                            float &pos, 
                            const Point pXtal) const {
  
  // create Volume Identifier for segments 0 & 11 of this crystal
  // volId info snagged from 
  // http://www.slac.stanford.edu/exp/glast/ground/software/geometry/docs/identifiers/geoId-RitzId.shtml

  const CalUtil::XtalIdx xtalIdx(xtalId);
  idents::VolumeIdentifier segm0Id, segm11Id;
  
  // init seg0 w/ info shared by both.
  segm0Id.append(m_eLATTowers);
  segm0Id.append(xtalIdx.getTwr().getRow());
  segm0Id.append(xtalIdx.getTwr().getCol());
  segm0Id.append(m_eTowerCAL);
  segm0Id.append(xtalIdx.getLyr().val());
  segm0Id.append(xtalIdx.getLyr().getDir().val()); 
  segm0Id.append(0);
  segm0Id.append(m_eXtal);

  segm11Id.append(m_eLATTowers);
  segm11Id.append(xtalIdx.getTwr().getRow());
  segm11Id.append(xtalIdx.getTwr().getCol());
  segm11Id.append(m_eTowerCAL);
  segm11Id.append(xtalIdx.getLyr().val());
  segm11Id.append(xtalIdx.getLyr().getDir().val()); 
 // segm11Id.append(xtalIdx.getCol().val());
  segm11Id.append(11);
  segm11Id.append(m_eXtal);

 
  // init segment specific info.
  segm0Id.append(5); // segment Id's set to middle of log
  segm11Id.append(5); 
  
  HepTransform3D transf;

  //get 3D transformation for segment 0 of this crystal
  m_detSvc->getTransform3DByID(segm0Id,&transf);
  //get position of the center of the segment 0
  const Vector vect0 = transf.getTranslation();

  //get 3D transformation for the last segment of this crystal
  m_detSvc->getTransform3DByID(segm11Id,&transf);
  //get position of the center of the last segment
  const Vector vect11 = transf.getTranslation();
          
  Point p0(0., 0., 0.);
  // position of the crystal center - need to bogus p0 to cast Vector -> Point
  const Point pCenter = p0 + (vect0+vect11)*0.5; 

  Point hit = p0 +(pXtal - pCenter);
  pos = (xtalId.isX()) ?  hit.y() : hit.x();
}
