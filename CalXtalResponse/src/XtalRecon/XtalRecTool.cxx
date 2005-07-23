// LOCAL
#include "XtalRecTool.h"

// GLAST
#include "idents/VolumeIdentifier.h"

// EXTLIB
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/MsgStream.h"
#include "CLHEP/Geometry/Transform3D.h"

// STD
#include <cmath>
#include <algorithm>

using namespace CalDefs;
using Event::CalDigi;
using Event::CalXtalRecData;

const char *XtalRecTool::m_tupleDesc = 
"RunID/i"
":EventID/i"
":adc[2]/s"
":adcPed[2]/F"
":ene/F"
":faceSignal[2]/F"
":asymCtr[2]/F"
":pos/F"
":dac[2]/F"
":asym/F"
":meanDAC/F"
":ped[2]/F"
":pedSig[2]/F"
":lacThresh[2]/F"
":mpd[2]/F"
":h1Limit[2]/F"
":twr/b"
":lyr/b"
":col/b"
":success/b"
":rng[2]/b"
":belowThresh[2]/b"
":xtalBelowThresh/b"
":saturated[2]/b"
":diode[2]/b"
;

static ToolFactory<XtalRecTool> s_factory;
const IToolFactory& XtalRecToolFactory = s_factory;

XtalRecTool::XtalRecTool( const string& type,
                          const string& name,
                          const IInterface* parent)
  : AlgTool(type,name,parent),
    m_calCalibSvc(0),
    m_tuple(0),
    m_tupleFile(0),
    m_tupleBranch(0),
    m_CsILength(0),
    m_eLATTowers(-1),
    m_eTowerCAL(-1),
    m_eXtal(-1),
    m_nCsISeg(-1)
{
  declareInterface<IXtalRecTool>(this);

  declareProperty("CalCalibSvc", m_calCalibSvcName="CalCalibSvc");
  declareProperty("tupleFilename",   m_tupleFilename="");
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

  // obtain CalCalibSvc
  sc = service(m_calCalibSvcName.value(), m_calCalibSvc);
  if (sc.isFailure()) {
    msglog << MSG::ERROR << "can't get CalCalibSvc." << endreq;
    return sc;
  }

  //-- EXTRACT INT CONSTANTS --//
  double value;
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
    if(!m_detSvc->getNumericConstByName((*iter).second, &value)) {
      MsgStream msglog(msgSvc(), name());
      msglog << MSG::ERROR << " constant " <<(*iter).second
             <<" not defined" << endreq;
      return StatusCode::FAILURE;
    } else *((*iter).first)= int(value); // store retrieved value 
  }

  //-- EXTRACT DOUBLE CONSTANTS --//

  // map containing pointers to double constants to be read
  // with their symbolic names from xml file used as a key 
  typedef map<double*,string> DPARAMAP;
  DPARAMAP dparam; 

  dparam[&m_CsILength]  = string("CsILength");
    
  for(DPARAMAP::iterator dIter=dparam.begin(); dIter!=dparam.end();dIter++){
    if(!m_detSvc->getNumericConstByName((*dIter).second,(*dIter).first)) {
      MsgStream msglog(msgSvc(), name());
      msglog << MSG::ERROR << " constant " <<(*dIter).second << " not defined" << endreq;
      return StatusCode::FAILURE;
    } 
  }

  // open optional tuple file
  if (m_tupleFilename.value() != "") {
    m_tupleFile = new TFile(m_tupleFilename.value().c_str(),"RECREATE","XtalRecTuple");
    if (!m_tupleFile) {
      msglog << MSG::ERROR << "Unable to create TTree object: " << m_tupleFilename << endreq;
      return StatusCode::FAILURE;
    }

    m_tuple = new TTree("XtalRecTuple","XtalRecTuple");
    if (!m_tuple) {
      msglog << MSG::ERROR << "Unable to create tuple" << endreq;
      return StatusCode::FAILURE;
    }

    m_tupleBranch = m_tuple->Branch("XtalRecTuple",
                                    &m_dat,
                                    m_tupleDesc);
    if (!m_tupleBranch) {
      msglog << MSG::ERROR << "Couldn't create tuple branch" << endreq;
      return StatusCode::FAILURE;
    }
  }

  return StatusCode::SUCCESS;
}

StatusCode XtalRecTool::calculate(const Event::EventHeader &evtHdr,
                                  const CalDigi &digi,
                                  CalXtalRecData &xtalRec,
                                  bool &belowThreshP,   
                                  bool &belowThreshN,   
                                  bool &xtalBelowThresh,
                                  bool &saturatedP,     
                                  bool &saturatedN,
                                  CalTupleEntry *calTupleEnt
                                  ) {
  StatusCode sc;

  // clear algorithm vars
  m_dat.Clear();
  
  // used for global access routines.
  CalXtalId xtalId = digi.getPackedId();
  m_dat.twr = xtalId.getTower();
  m_dat.lyr = xtalId.getLayer();
  m_dat.col = xtalId.getColumn();

  // currently allways using 1st readout
  CalDigi::CalXtalReadoutCol::const_iterator ro = digi.getReadoutCol().begin();

  // get readout range number for both crystal faces
  m_dat.rng[POS_FACE] = (*ro).getRange(CalXtalId::POS); 
  m_dat.rng[NEG_FACE] = (*ro).getRange(CalXtalId::NEG); 
  
  // get adc values 
  m_dat.adc[POS_FACE] = (*ro).getAdc(CalXtalId::POS);   
  m_dat.adc[NEG_FACE] = (*ro).getAdc(CalXtalId::NEG); 

  // initialize return values
  xtalBelowThresh = 
    belowThreshP = belowThreshN = 
    saturatedP = saturatedN =
    false;

  ///////////////////////////////////////
  //-- STEP 0: RETRIEVE CALIBRATIONS --//
  ///////////////////////////////////////
  sc = retrieveCalib(xtalId);
  if (sc.isFailure()) return sc;
 
  //////////////////////////////////////
  //-- STEP 1: PEDESTAL SUBTRACTION --//
  //////////////////////////////////////
  m_dat.adcPed[POS_FACE] = m_dat.adc[POS_FACE] - m_dat.ped[POS_FACE];   // ped subtracted ADC
  m_dat.adcPed[NEG_FACE] = m_dat.adc[NEG_FACE] - m_dat.ped[NEG_FACE]; // ped subtracted ADC

  
  /////////////////////////////////
  //-- STEP 2: NOISE REDUCTION --//
  /////////////////////////////////
  
  // LEX8 range is compared against 0.5 * lac threshold
  // we throw out entire xtal if adc is too low 
  if (m_dat.rng[POS_FACE] == LEX8 && m_dat.adcPed[POS_FACE] < m_dat.lacThresh[POS_FACE]*0.5)
    belowThreshP = true, xtalBelowThresh = true;
  if (m_dat.rng[NEG_FACE] == LEX8 && m_dat.adcPed[NEG_FACE] < m_dat.lacThresh[NEG_FACE]*0.5)
    belowThreshN = true, xtalBelowThresh = true;

  // other ranges are compared against 5 sigma energy 
  if (m_dat.rng[POS_FACE] != LEX8 && m_dat.adcPed[POS_FACE] < m_dat.pedSig[POS_FACE]*5.0)
    belowThreshP = true;
  if (m_dat.rng[NEG_FACE] != LEX8 && m_dat.adcPed[NEG_FACE] < m_dat.pedSig[NEG_FACE]*5.0)
    belowThreshN = true;
  

  ///////////////////////////////////////////
  //-- STEP 3: CONVERT ADCs -> DAC SCALE --//
  ///////////////////////////////////////////

  for (FaceNum face; face.isValid(); face++) {
    double tmp;
    sc = m_calCalibSvc->evalDAC(CalXtalId(m_dat.twr,m_dat.lyr,m_dat.col,face,m_dat.rng[face]), 
                                m_dat.adcPed[face], tmp);
    m_dat.dac[face] = tmp;
    if (sc.isFailure()) return sc;
  }
  
  // check for invalid dac vals (i need to take the sqrt AND the log)
  if (m_dat.dac[POS_FACE] <= 0 || m_dat.dac[NEG_FACE] <= 0) {
    // create MsgStream only when needed for performance
    MsgStream msglog(msgSvc(), name()); 
    msglog << MSG::VERBOSE;
    // need to use .stream() to get xtalId to pretty print.
    msglog.stream() << "minimal signal... DAC val <= 0, can't calculate energy/position." 
                    <<  " xtal=[" << xtalId << ']';
    msglog << endreq;
    return StatusCode::FAILURE;
  }

  m_dat.diode[POS_FACE] = RngNum(m_dat.rng[POS_FACE]).getDiode();
  m_dat.diode[NEG_FACE] = RngNum(m_dat.rng[NEG_FACE]).getDiode();

   

  ////////////////////////////////////
  //-- STEP 4: CALCULATE POSITION --//
  ////////////////////////////////////

  m_dat.asym = log(m_dat.dac[POS_FACE]/m_dat.dac[NEG_FACE]);
  double tmp;
  if (m_dat.diode[POS_FACE] == LRG_DIODE)
    if (m_dat.diode[NEG_FACE] == LRG_DIODE)
      sc = m_calCalibSvc->evalPosLrg(xtalId, m_dat.asym, tmp);   
    else
      sc = m_calCalibSvc->evalPosNSPB(xtalId, m_dat.asym, tmp);
  else
    if (m_dat.diode[NEG_FACE] == SM_DIODE)
      sc = m_calCalibSvc->evalPosSm(xtalId, m_dat.asym, tmp);
    else
      sc = m_calCalibSvc->evalPosPSNB(xtalId, m_dat.asym, tmp);
  m_dat.pos = tmp;
  if (sc.isFailure()) return sc;

  // generate position vector from scalar longitudinal position
  Point pXtal;
  pos2Point(m_dat.pos, pXtal);


  ///////////////////////////////////////////
  //-- STEP 5: CALCULATE XTAL-FACE SIGNAL--//
  ///////////////////////////////////////////

  // 1st multiply each dac val by overall gain.
  m_dat.faceSignal[POS_FACE] = m_dat.dac[POS_FACE]*m_dat.mpd[m_dat.diode[POS_FACE]];
  m_dat.faceSignal[NEG_FACE] = m_dat.dac[NEG_FACE]*m_dat.mpd[m_dat.diode[NEG_FACE]];
         
  // 2nd correct for overall asymmetry of diodes (asym at center of xtal)       
  m_dat.faceSignal[POS_FACE] *= exp(-1*m_dat.asymCtr[m_dat.diode[POS_FACE]]/2);
  m_dat.faceSignal[NEG_FACE] *= exp(   m_dat.asymCtr[m_dat.diode[NEG_FACE]]/2);


  ///////////////////////////////////////
  //-- STEP 6: CALCULATE TRUE ENERGY --//
  ///////////////////////////////////////
  
  //-- convert diodes to same size (if needed)
  DiodeNum mpdDiode = m_dat.diode[POS_FACE]; // diode size to use for MeVPerDAC conversion
  CalVec<DiodeNum,float> mpdDac(DiodeNum::N_VALS);
  mpdDac[POS_FACE]  = m_dat.dac[POS_FACE];   // dac value to use for MeVPerDAC conversion
  mpdDac[NEG_FACE]  = m_dat.dac[NEG_FACE];   // dac value to use for MeVPerDAC conversion

  // if diodes on each face differ convert to small on both faces
  if (m_dat.diode[POS_FACE] != m_dat.diode[NEG_FACE]) {
    if (m_dat.diode[POS_FACE] == LRG_DIODE)
      sc = largeDAC2Small(POS_FACE,mpdDac[POS_FACE],mpdDac[POS_FACE]);
    else // m_dat.diode[NEG_FACE] == LRG_DIODE
      sc = largeDAC2Small(NEG_FACE,mpdDac[NEG_FACE],mpdDac[NEG_FACE]);
    
    if (sc.isFailure()) return sc;
    
    mpdDiode = SM_DIODE;
  }

  m_dat.meanDAC = sqrt(mpdDac[POS_FACE]*mpdDac[NEG_FACE]);
  m_dat.ene = m_dat.meanDAC*m_dat.mpd[mpdDiode];



  /////////////////////////////////////////////
  //-- STEP 7: DETECT SATURATED HEX1 RANGE --//
  /////////////////////////////////////////////
  if (m_dat.rng[POS_FACE] == HEX1)
    if (m_dat.adc[POS_FACE] >= m_dat.h1Limit[POS_FACE]) saturatedP = true;
  if (m_dat.rng[NEG_FACE] == HEX1)
    if (m_dat.adc[NEG_FACE] >= m_dat.h1Limit[NEG_FACE]) saturatedN = true;



  ////////////////////////////////////
  //-- STEP 8: POPULATE TDS CLASS --//
  ////////////////////////////////////
  CalXtalRecData::CalRangeRecData rngRec(m_dat.rng[POS_FACE], m_dat.ene, 
                                         m_dat.rng[NEG_FACE], m_dat.ene);
  rngRec.setPosition(pXtal);

  xtalRec.addRangeRecData(rngRec);
  
  
  //////////////////////////////////////////////////
  //-- STEP 9: POPULATE XtalRecTuple (OPTIONAL) --//
  //////////////////////////////////////////////////

  if (m_tuple) {
    m_dat.RunID                   = evtHdr.run();
    m_dat.EventID                   = evtHdr.event();
    m_dat.success                = true;
    m_dat.belowThresh[POS_FACE]  = belowThreshP;
    m_dat.belowThresh[NEG_FACE]  = belowThreshN;
    m_dat.saturated[POS_FACE]    = saturatedP;
    m_dat.saturated[NEG_FACE]    = saturatedN;
    m_dat.xtalBelowThresh        = xtalBelowThresh;

    m_tuple->Fill();
  }


  ///////////////////////////////////////////////
  //-- STEP 10: POPULATE CalTuple (OPTIONAL) --//
  ///////////////////////////////////////////////
  if (calTupleEnt) {
    calTupleEnt->m_calXtalAdcPed[m_dat.twr][m_dat.lyr][m_dat.col][POS_FACE] = m_dat.adcPed[POS_FACE];
    calTupleEnt->m_calXtalAdcPed[m_dat.twr][m_dat.lyr][m_dat.col][NEG_FACE] = m_dat.adcPed[NEG_FACE];
    
    calTupleEnt->m_calXtalFaceSignal[m_dat.twr][m_dat.lyr][m_dat.col][POS_FACE] = m_dat.faceSignal[POS_FACE];
    calTupleEnt->m_calXtalFaceSignal[m_dat.twr][m_dat.lyr][m_dat.col][NEG_FACE] = m_dat.faceSignal[NEG_FACE];
  }

  return StatusCode::SUCCESS;
}

StatusCode XtalRecTool::largeDAC2Small(FaceNum face, float largeDAC, float &smallDAC) {
  StatusCode sc;

  if (face == POS_FACE) {
    // STRATEGY:
    // convert PL -> PS
    // AsymNSPB = log(PL / NS)
    // AsymSm   = log(PS / NS)
    // exp(asymSm)/exp(asymNSPB) = (PS/NS)/(PL/NS) = PS/PL
    //      exp(asymSm-asymNSPB) = PS/PL
    //                        PS = PL*exp(asymSm=asymNSPB
      
    double asymNSPB;    
    sc = m_calCalibSvc->evalAsymNSPB(CalXtalId(m_dat.twr,m_dat.lyr,m_dat.col), 
                                     0.0, asymNSPB);
    if (sc.isFailure()) return sc;

    smallDAC = largeDAC * exp(m_dat.asymCtr[SM_DIODE] - asymNSPB);
  } 
  else { 
    // STRATEGY:
    // convert NL -> NS
    // asymPSNB = log(PS/NL)
    // asymSm   = log(PS/NS)
    // exp(asymPSNB)/exp(asymSm) = (PS/NL)/(PS/NS)
    //      exp(asymPSNB-asymSm) = NS/NL
    //                        NS = NL*exp(asymPSNB-asymSm)

    double asymPSNB;
    sc = m_calCalibSvc->evalAsymPSNB(CalXtalId(m_dat.twr,m_dat.lyr,m_dat.col), 
                                     0.0, asymPSNB);
    if (sc.isFailure()) return sc;

    smallDAC = largeDAC * exp(asymPSNB - m_dat.asymCtr[SM_DIODE]);
  }

  return StatusCode::SUCCESS;

}


/** \brief convert longitudinal pos from XtalReconTool to 3d point in LAT geometry space

\param pXtal ouput position vector
\param pos input longitudinal position in mm from center of xtal
*/
void XtalRecTool::pos2Point(float pos, Point &pXtal) {
  //-- CONSTRUCT GEOMETRY VECTORS FOR XTAL --//

  TwrNum twr(m_dat.twr);
  LyrNum lyr(m_dat.lyr);

  // create Volume Identifier for segment 0 of this crystal
  // volId info snagged from 
  // http://www.slac.stanford.edu/exp/glast/ground/software/geometry/docs/identifiers/geoId-RitzId.shtml
  idents::VolumeIdentifier segm0Id;
  segm0Id.append(m_eLATTowers);
  segm0Id.append(twr.getRow());
  segm0Id.append(twr.getCol());
  segm0Id.append(m_eTowerCAL);
  segm0Id.append(lyr);
  segm0Id.append(lyr.getDir()); 
  segm0Id.append(m_dat.col);
  segm0Id.append(m_eXtal);
  segm0Id.append(0); // segment Id

  HepTransform3D transf;

  //get 3D transformation for segment 0 of this crystal
  m_detSvc->getTransform3DByID(segm0Id,&transf);
  //get position of the center of the segment 0
  Vector vect0 = transf.getTranslation();
  // create Volume Identifier for the last segment of this crystal
  idents::VolumeIdentifier segm11Id;
  // copy all fields from segm0Id, except segment number
  for(int ifield = 0; ifield < fSegment; ifield++)
    segm11Id.append(segm0Id[ifield]);
  segm11Id.append(m_nCsISeg-1); // set segment number for the last segment
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

StatusCode XtalRecTool::retrieveCalib(const CalXtalId &xtalId) {
  StatusCode sc;

  //-- RETRIEVE MEV PER DAC--// 
  CalibData::ValSig mpdLrg, mpdSm;
  sc = m_calCalibSvc->getMeVPerDac(xtalId, mpdLrg, mpdSm);
  if (sc.isFailure()) return sc;
  m_dat.mpd[LRG_DIODE] = mpdLrg.getVal();
  m_dat.mpd[SM_DIODE]  = mpdSm.getVal();

  //-- RETRIEVE ASYM @ CTR OF XTAL --//
  double tmp;
  sc = m_calCalibSvc->evalAsymLrg(xtalId, 0, tmp);
  m_dat.asymCtr[LRG_DIODE] = tmp;
  if (sc.isFailure()) return sc;
  
  sc = m_calCalibSvc->evalAsymSm(xtalId, 0, tmp);
  m_dat.asymCtr[SM_DIODE] = tmp;
  if (sc.isFailure()) return sc;

  for (FaceNum face; face.isValid(); face++) {
    // pedestals
    float tmp;
    CalXtalId rngXtalId(m_dat.twr,m_dat.lyr,m_dat.col, face, m_dat.rng[face]);
    sc = m_calCalibSvc->getPed(rngXtalId, m_dat.ped[face], m_dat.pedSig[face], tmp);
    if (sc.isFailure()) return sc;

    // Threshold constants
    CalibData::ValSig fle,fhe,lac;
    CalXtalId faceXtalId(m_dat.twr, m_dat.lyr, m_dat.col, face);
    sc = m_calCalibSvc->getTholdCI(faceXtalId,fle,fhe,lac);
    m_dat.lacThresh[face] = lac.getVal();

    //-- RETRIEVE HEX1 ULD --//
    CalibData::ValSig uldThold;
    CalXtalId h1XtalId(m_dat.twr,m_dat.lyr,m_dat.col,face,HEX1);
    sc = m_calCalibSvc->getULDCI(h1XtalId,uldThold);
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
