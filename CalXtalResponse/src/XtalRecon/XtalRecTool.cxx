// $Header$
/** @file
    @author Zach Fewtrell
 */

// LOCAL
#include "XtalRecTool.h"

// GLAST
#include "idents/VolumeIdentifier.h"
#include "CalUtil/CalArray.h"

// EXTLIB
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/MsgStream.h"
#include "CLHEP/Geometry/Transform3D.h"

// STD
#include <cmath>
#include <map>
#include <algorithm>

using namespace CalUtil;
using Event::CalDigi;
using Event::CalXtalRecData;

static ToolFactory<XtalRecTool> s_factory;
const IToolFactory& XtalRecToolFactory = s_factory;

XtalRecTool::XtalRecTool( const string& type,
                          const string& name,
                          const IInterface* parent)
  : AlgTool(type,name,parent),
    m_calCalibSvc(0),
    m_tuple(0),
    m_tupleFile(0),
    m_CsILength(0),
    m_eLATTowers(-1),
    m_eTowerCAL(-1),
    m_eXtal(-1),
    m_nCsISeg(-1)
{
  declareInterface<IXtalRecTool>(this);

  declareProperty("CalCalibSvc", m_calCalibSvcName="CalCalibSvc");
  declareProperty("tupleFilename",  m_tupleFilename="");
}

StatusCode XtalRecTool::initialize() {
  MsgStream msglog(msgSvc(), name());
  msglog << MSG::INFO << "initialize" << endreq;

  StatusCode sc;

  //-- jobOptions --//
  if ((sc = setProperties()).isFailure()) {
    msglog << MSG::ERROR << "Failed to set properties" << endreq;
    return sc;
  }

  //-- tuple --//
  if (m_tupleFilename.value().length() > 0) {
    m_tupleFile = new TFile(m_tupleFilename.value().c_str(),"RECREATE","XtalRecTuple");
    if (!m_tupleFile)
      // allow to continue w/out tuple file as it is not a _real_ failure
      msglog << MSG::ERROR << "Unable to create TFile object: " << m_tupleFilename << endreq;
    else {
      m_tuple = new TTree("XtalRecTuple","XtalRecTuple");
      if (!m_tuple) {
        msglog << MSG::ERROR << "Unable to create tuple" << endreq;
        return StatusCode::FAILURE;
      }

      if (!m_tuple->Branch("RunID",           &m_dat.RunID,              "RunID/i")           ||
          !m_tuple->Branch("EventID",         &m_dat.EventID,            "EventID/i")         ||
          !m_tuple->Branch("adc",             m_dat.adc.begin(),         "adc[2]/s")          ||
          !m_tuple->Branch("adcPed",          m_dat.adcPed.begin(),      "adcPed[2]/F")       ||
          !m_tuple->Branch("ene",             &m_dat.ene,                "ene/F")             ||
          !m_tuple->Branch("faceSignal",      m_dat.faceSignal.begin(),  "faceSignal[2]/F")   ||
          !m_tuple->Branch("asymCtr",         m_dat.asymCtr.begin(),     "asymCtr[2]/F")      ||
          !m_tuple->Branch("pos",             &m_dat.pos,                "pos/F")             ||
          !m_tuple->Branch("dac",             m_dat.dac.begin(),         "dac[2]/F")          ||
          !m_tuple->Branch("asym",            &m_dat.asym,               "asym/F")            ||
          !m_tuple->Branch("meanDAC",         &m_dat.meanDAC,            "meanDAC/F")         ||
          !m_tuple->Branch("ped",             m_dat.ped.begin(),         "ped[2]/F")          ||
          !m_tuple->Branch("pedSig",          m_dat.pedSig.begin(),      "pedSig[2]/F")       ||
          !m_tuple->Branch("lacThresh",       m_dat.lacThresh.begin(),   "lacThresh[2]/F")    ||
          !m_tuple->Branch("mpd",             m_dat.mpd.begin(),         "mpd[2]/F")          ||
          !m_tuple->Branch("h1Limit",         m_dat.h1Limit.begin(),     "h1Limit[2]/F")      ||
          !m_tuple->Branch("twr",             &m_dat.twr,                "twr/b")             ||
          !m_tuple->Branch("lyr",             &m_dat.lyr,                "lyr/b")             ||
          !m_tuple->Branch("col",             &m_dat.col,                "col/b")             ||
          !m_tuple->Branch("rng",             m_dat.rng.begin(),         "rng[2]/b")          ||
          !m_tuple->Branch("belowThresh",     (void*)m_dat.belowThresh.begin(), "belowThresh[2]/b")  ||
          !m_tuple->Branch("xtalBelowThresh", (void*)&m_dat.xtalBelowThresh,    "xtalBelowThresh/b") ||
          !m_tuple->Branch("saturated",       (void*)m_dat.saturated.begin(),   "saturated[2]/b")    ||
          !m_tuple->Branch("diode",           m_dat.diode.begin(),       "diode[2]/b")) {
        msglog << MSG::ERROR << "Couldn't create tuple branch" << endreq;
        return StatusCode::FAILURE;
      } 
    }
    
  } // optional tuple
  
  // obtain CalCalibSvc
  sc = service(m_calCalibSvcName.value(), m_calCalibSvc);
  if (sc.isFailure()) {
    msglog << MSG::ERROR << "can't get CalCalibSvc." << endreq;
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
  m_CsILength = tmp;
  if (sc.isFailure()) {
    msglog << MSG::ERROR << " constant CsILength not defined" << endreq;
    return StatusCode::FAILURE;
  }

  return StatusCode::SUCCESS;
}

/**
   - convert adc readouts to DAC scale
   - compute centroid position & intesity from DAC scale
   - check for noise threshold & xtal saturation
*/
StatusCode XtalRecTool::calculate(const Event::CalDigi &digi,
                                  Event::CalXtalRecData &xtalRec,
                                  CalArray<FaceNum, bool> &belowThresh,
                                  bool &xtalBelowThresh,
                                  CalArray<FaceNum, bool> &saturated,
                                  const Event::EventHeader *evtHdr) {
  StatusCode sc;

  // clear algorithm vars
  m_dat.Clear();

  // initialize return values
  belowThresh.fill(false);
  saturated.fill(false);
  xtalBelowThresh = false;
  
  // used for global access routines.
  m_dat.xtalIdx = XtalIdx(digi.getPackedId());
  m_dat.twr = m_dat.xtalIdx.getTwr();
  m_dat.lyr = m_dat.xtalIdx.getLyr();
  m_dat.col = m_dat.xtalIdx.getCol();



  // currently allways using 1st readout
  CalDigi::CalXtalReadoutCol::const_iterator ro = 
    digi.getReadoutCol().begin();

  // get readout range number for both crystal faces
  for (FaceNum face; face.isValid();  face++) {
    // get adc range
    m_dat.rng[face] = (*ro).getRange(face); 
    // get adc values 
    m_dat.adc[face] = (*ro).getAdc(face);   
  }

  ///////////////////////////////////////
  //-- STEP 0: RETRIEVE CALIBRATIONS --//
  ///////////////////////////////////////
  sc = retrieveCalib();
  if (sc.isFailure()) return sc;
 
  //////////////////////////////////////
  //-- STEP 1: PEDESTAL SUBTRACTION --//
  //////////////////////////////////////
  for (FaceNum face; face.isValid();  face++) {
    m_dat.adcPed[face] = 
      m_dat.adc[face] - m_dat.ped[face]; // ped subtracted ADC

  
    /////////////////////////////////
    //-- STEP 2: NOISE REDUCTION --//
    /////////////////////////////////
    
    // LEX8 range is compared against 0.5 * lac threshold
    // we throw out entire xtal if adc is too low 
    if (m_dat.rng[face] == LEX8 && m_dat.adcPed[face] < m_dat.lacThresh[face]*0.5)
      belowThresh[face] = true, xtalBelowThresh = true;

    // other ranges are compared against 5 sigma energy 
    if (m_dat.rng[face] != LEX8 && m_dat.adcPed[face] < m_dat.pedSig[face]*5.0)
      belowThresh[face] = true;

    ///////////////////////////////////////////
    //-- STEP 3: CONVERT ADCs -> DAC SCALE --//
    ///////////////////////////////////////////
	RngIdx rngIdx(m_dat.xtalIdx, face, m_dat.rng[face]);
    sc = m_calCalibSvc->evalDAC(rngIdx, m_dat.adcPed[face], m_dat.dac[face]);
    if (sc.isFailure()) return sc;
        
    // check for invalid dac vals (i need to take the sqrt AND the log)
    if (m_dat.dac[face] <= 0) {
      // create MsgStream only when needed (for performance)
      MsgStream msglog(msgSvc(), name()); 
      msglog << MSG::VERBOSE;
      // need to use .stream() to get xtalId to pretty print.
      msglog.stream() << "minimal signal... DAC val <= 0, can't calculate energy/position." 
                      <<  " xtal=[" << m_dat.xtalIdx.getCalXtalId() << ']';
      msglog << endreq;

      xtalBelowThresh = true;
      belowThresh[face] = true;

      // this is not a FAILURE condition, just not enough signal
      // return w/out populating xtalRec structure
      return StatusCode::SUCCESS;
    }
        


    m_dat.diode[face] = m_dat.rng[face].getDiode();

    /////////////////////////////////////////////
    //-- STEP 4: DETECT SATURATED HEX1 RANGE --//
    /////////////////////////////////////////////
    if (m_dat.rng[face] == HEX1)
      if (m_dat.adc[face] >= m_dat.h1Limit[face]) saturated[face] = true;
    
  
    ////////////////////////////////////////////
    //-- STEP 5: POPULATE TUPLES (OPTIONAL) --//
    ////////////////////////////////////////////
    
    // need face signal for either tuple
    if (m_tuple) {
      sc = m_calCalibSvc->evalFaceSignal(rngIdx,
                                         m_dat.adcPed[face], 
                                         m_dat.faceSignal[face]);
      if (sc.isFailure()) return sc;
    }
  }



  ////////////////////////////////////
  //-- STEP 7: CALCULATE POSITION --//
  ////////////////////////////////////

  m_dat.asym = log(m_dat.dac[POS_FACE]/m_dat.dac[NEG_FACE]);
  if (m_dat.diode[POS_FACE] == LRG_DIODE)
    if (m_dat.diode[NEG_FACE] == LRG_DIODE)
      sc = m_calCalibSvc->evalPosLrg(m_dat.xtalIdx, m_dat.asym, m_dat.pos);   
    else
      sc = m_calCalibSvc->evalPosNSPB(m_dat.xtalIdx, m_dat.asym, m_dat.pos);
  else
    if (m_dat.diode[NEG_FACE] == SM_DIODE)
      sc = m_calCalibSvc->evalPosSm(m_dat.xtalIdx, m_dat.asym, m_dat.pos);
    else
      sc = m_calCalibSvc->evalPosPSNB(m_dat.xtalIdx, m_dat.asym, m_dat.pos);
  if (sc.isFailure()) return sc;

  // 1st clip position to end of xtal
  m_dat.pos = min<float>(   m_CsILength/2,m_dat.pos);
  m_dat.pos = max<float>(-1*m_CsILength/2,m_dat.pos);
  
  // generate position vector from scalar longitudinal position
  Point pXtal;
  pos2Point(m_dat.pos, pXtal);


  ///////////////////////////////////////
  //-- STEP 8: CALCULATE TRUE ENERGY --//
  ///////////////////////////////////////
  
  //-- convert diodes to same size (if needed)
  DiodeNum mpdDiode(m_dat.diode[POS_FACE]); // diode size to use for MeVPerDAC conversion
  CalArray<FaceNum,float> mpdDac;
  mpdDac[POS_FACE]  = m_dat.dac[POS_FACE];   // dac value to use for MeVPerDAC conversion
  mpdDac[NEG_FACE]  = m_dat.dac[NEG_FACE];   // dac value to use for MeVPerDAC conversion

  // if diodes on each face differ convert to small on both faces
  if (m_dat.diode[POS_FACE].getInt() != m_dat.diode[NEG_FACE].getInt()) {
    if (m_dat.diode[POS_FACE] == LRG_DIODE)
      sc = largeDAC2Small(POS_FACE, m_dat.pos, mpdDac[POS_FACE], mpdDac[POS_FACE]);
    else // m_dat.diode[NEG_FACE] == LRG_DIODE
      sc = largeDAC2Small(NEG_FACE, m_dat.pos, mpdDac[NEG_FACE], mpdDac[NEG_FACE]);
    
    if (sc.isFailure()) return sc;
    
    mpdDiode = SM_DIODE;
  }

  m_dat.meanDAC = sqrt(mpdDac[POS_FACE]*mpdDac[NEG_FACE]);
  m_dat.ene = m_dat.meanDAC*m_dat.mpd[mpdDiode];


  /////////////////////////////////////////////////
  //-- STEP 9: POPULATE OPTIONAL XTALRECTUPLE  --//
  /////////////////////////////////////////////////
  if (m_tuple) {
    if (evtHdr) {
      m_dat.RunID                   = evtHdr->run();
      m_dat.EventID                 = evtHdr->event();
    }
    copy(belowThresh.begin(), belowThresh.end(), m_dat.belowThresh.begin());
    copy(saturated.begin(), saturated.end(),   m_dat.saturated.begin());
    m_dat.xtalBelowThresh        = xtalBelowThresh;

    //-- RETRIEVE ASYM @ CTR OF XTAL --//
    sc = m_calCalibSvc->evalAsymLrg(m_dat.xtalIdx, 0, m_dat.asymCtr[LRG_DIODE]);
    if (sc.isFailure()) return sc;
  
    sc = m_calCalibSvc->evalAsymSm(m_dat.xtalIdx, 0, m_dat.asymCtr[SM_DIODE]);
    if (sc.isFailure()) return sc;

    m_tuple->Fill();
  }


  ////////////////////////////////////
  //-- STEP 10: POPULATE TDS CLASS --//
  ////////////////////////////////////
  CalXtalRecData::CalRangeRecData rngRec((CalXtalId::AdcRange)m_dat.rng[POS_FACE], m_dat.ene, 
                                         (CalXtalId::AdcRange)m_dat.rng[NEG_FACE], m_dat.ene);
  rngRec.setPosition(pXtal);

  xtalRec.addRangeRecData(rngRec);

  return StatusCode::SUCCESS;
}

StatusCode XtalRecTool::largeDAC2Small(FaceNum face, float pos, float largeDAC, float &smallDAC) {
  StatusCode sc;

  // small diode asymmetry is used for both if() cases
  float asymSm;
  sc = m_calCalibSvc->evalAsymSm(m_dat.xtalIdx, pos, asymSm);
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
    sc = m_calCalibSvc->evalAsymNSPB(m_dat.xtalIdx, pos, asymNSPB);
    if (sc.isFailure()) return sc;

    smallDAC = largeDAC * exp(asymSm - asymNSPB);
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
    sc = m_calCalibSvc->evalAsymPSNB(m_dat.xtalIdx, pos, asymPSNB);
    if (sc.isFailure()) return sc;

    smallDAC = largeDAC * exp(asymPSNB - asymSm);
  }

  return StatusCode::SUCCESS;

}


/** \brief convert longitudinal pos along cal xtal to 3d point in LAT geometry space

\param pXtal ouput position vector
\param pos input longitudinal position in mm from center of xtal
*/
void XtalRecTool::pos2Point(float pos, Point &pXtal) {
  //-- CONSTRUCT GEOMETRY VECTORS FOR XTAL --//

  // create Volume Identifier for segments 0 & 11 of this crystal
  
  // volId info snagged from 
  // http://www.slac.stanford.edu/exp/glast/ground/software/geometry/docs/identifiers/geoId-RitzId.shtml
  idents::VolumeIdentifier segm0Id, segm11Id;
  
  // init seg0 w/ info shared by both.
  segm0Id.append(m_eLATTowers);
  segm0Id.append(m_dat.twr.getRow());
  segm0Id.append(m_dat.twr.getCol());
  segm0Id.append(m_eTowerCAL);
  segm0Id.append(m_dat.lyr);
  segm0Id.append(m_dat.lyr.getDir()); 
  segm0Id.append(m_dat.col);
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
  // 'pos' is in units of xtal Length, convert to rel units (-1->1)
  pos /= m_CsILength; 
  pXtal = pCenter+dirXtal*pos;
}

StatusCode XtalRecTool::retrieveCalib() {
  StatusCode sc;

  //-- RETRIEVE MEV PER DAC--// 
  sc = m_calCalibSvc->getMPD(m_dat.xtalIdx, m_dat.mpd);
  if (sc.isFailure()) return sc;

  for (FaceNum face; face.isValid(); face++) {
    FaceIdx faceIdx(m_dat.twr, m_dat.lyr, m_dat.col, face);
    RngIdx rngIdx(m_dat.twr, 
                  m_dat.lyr,
                  m_dat.col, 
                  face, m_dat.rng[face]);

    // pedestals
    float cos; // not used
    sc = m_calCalibSvc->getPed(rngIdx, m_dat.ped[face], m_dat.pedSig[face], cos);
    if (sc.isFailure()) return sc;

    // Threshold constants
    CalibData::ValSig fle, fhe, lac;
    sc = m_calCalibSvc->getTholdCI(faceIdx,fle,fhe,lac);
    m_dat.lacThresh[face] = lac.getVal();

    //-- RETRIEVE HEX1 ULD --//
    CalibData::ValSig uldThold;
    sc = m_calCalibSvc->getULDCI(RngIdx(faceIdx, HEX1),uldThold);
    if (sc.isFailure()) return sc;
    m_dat.h1Limit[face] = uldThold.getVal();
  }  

  return StatusCode::SUCCESS;
}

StatusCode XtalRecTool::finalize() {
  // make sure optional tuple is closed out                                                        
  if (m_tupleFile) {
    m_tupleFile->Write();
    m_tupleFile->Close(); // trees deleted                                                         
    delete m_tupleFile;
    m_tupleFile = 0;
  }

  return StatusCode::SUCCESS;
}
