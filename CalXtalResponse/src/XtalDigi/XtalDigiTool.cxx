// Include files

// LOCAL
#include "XtalDigiTool.h"

// GLAST
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"

// EXTLIB
#include "GaudiKernel/ToolFactory.h"

#include "GaudiKernel/MsgStream.h"

#include "CLHEP/Random/RandGauss.h"

// STD

using namespace CalDefs;


static ToolFactory<XtalDigiTool> s_factory;
const IToolFactory& XtalDigiToolFactory = s_factory;

const char *XtalDigiTool::m_tupleDesc = 
"RunID/i"
":EventID/i"
":adc[2][4]/F"
":adcPed[2][4]/F"
":nMCHits/i"
":nCsIHits/i"
":nDiodeHits[2][2]/i"
":sumEneCsI/F"
":sumEne/F"
":diodeDAC[2][2]/F"
":csiWeightedPos/F"
":ped[2][4]/F"
":pedSig[2][4]/F"
":lacThresh[2]/F"
":fleThresh[2]/F"
":fheThresh[2]/F"
":uldThold[2][4]/F"
":mpd[2]/F"
":twr/b"
":lyr/b"
":col/b"
":success/b"
":rng[2]/b"
":lac[2]/b"
":fle[2]/b"
":fhe[2]/b"
":saturated[2]/b"
;


XtalDigiTool::XtalDigiTool( const string& type,
                            const string& name,
                            const IInterface* parent)
  : AlgTool(type,name,parent),
    m_calCalibSvc(0),
    m_FailSvc(0),
    m_CsILength(0),
    m_eDiodeMSm(0),
    m_eDiodePSm(0),
    m_eDiodeMLarge(0),
    m_eDiodePLarge(0),
    m_ePerMeVInDiode(0),
    m_maxAdc(-1),
    m_tuple(0),
    m_tupleFile(0),
    m_tupleBranch(0)
{
  declareInterface<IXtalDigiTool>(this);

  declareProperty("CalCalibSvc",        m_calCalibSvcName         = "CalCalibSvc");
  declareProperty("NoRandomNoise",      m_noRandNoise             = false);
  declareProperty("tupleFilename",          m_tupleFilename               = "");
  declareProperty("AllowNoiseOnlyTrig", m_allowNoiseOnlyTrig = true);
}

StatusCode XtalDigiTool::initialize() {
  MsgStream msglog(msgSvc(), name());   
  msglog << MSG::INFO << "initialize" << endreq;

  StatusCode sc;

  double val;
  typedef map<int*,string> PARAMAP;

  //-- jobOptions --//
  if ((sc = setProperties()).isFailure()) {
    msglog << MSG::ERROR << "Failed to set properties" << endreq;
    return sc;
  }

  // retrieve cal failure service
  sc = service("CalFailureModeSvc", m_FailSvc);
  if (sc.isFailure() ) {
    msglog << MSG::INFO << "  Did not find CalFailureMode service" << endreq;
    m_FailSvc = 0;
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
  param[&m_nCsISeg]      = string("nCsISeg");
  param[&m_eXtal]        = string("eXtal");
  param[&m_eDiodeMSm]    = string("eDiodeMSmall");
  param[&m_eDiodePSm]    = string("eDiodePSmall");
  param[&m_eDiodeMLarge] = string("eDiodeMLarge");
  param[&m_eDiodePLarge] = string("eDiodePLarge");
  param[m_ePerMeV+1]     = string("cal.ePerMeVSmall");
  param[m_ePerMeV]       = string("cal.ePerMevLarge");
  param[&m_maxAdc]       = string("cal.maxAdcValue");

  for(PARAMAP::iterator it=param.begin(); it!=param.end();it++)
    if(!detSvc->getNumericConstByName((*it).second, &val)) {
      msglog << MSG::ERROR << " constant " <<(*it).second <<" not defined" << endreq;
      return StatusCode::FAILURE;
    } else *((*it).first)=(int)val;

  // doubles are done separately
  double tmp;
  sc = detSvc->getNumericConstByName("CsILength", &tmp);
  m_CsILength = tmp;
  if (sc.isFailure()) {
    msglog << MSG::ERROR << " constant CsILength not defined" << endreq;
    return StatusCode::FAILURE;
  }
  
  // eperMeVInDiode originally from CalDigi XML file
  // but i don't want dependency, so i'm hard coding it.  (just a testing tool anyway :)
  m_ePerMeVInDiode = 2.77e5;
  
  // obtain CalCalibSvc
  sc = service(m_calCalibSvcName.value(), m_calCalibSvc);
  if (sc.isFailure()) {
    msglog << MSG::ERROR << "can't get CalCalibSvc." << endreq;
    return sc;
  }

  // open optional tuple file
  if (m_tupleFilename.value() != "") {
    m_tupleFile = new TFile(m_tupleFilename.value().c_str(),"RECREATE","XtalDigiTuple");
    if (!m_tupleFile) {
      msglog << MSG::ERROR << "Unable to create TTree object: " << m_tupleFilename << endreq;
      return StatusCode::FAILURE;
    }

    m_tuple = new TTree("XtalDigiTuple","XtalDigiTuple");
    if (!m_tuple) {
      msglog << MSG::ERROR << "Unable to create tuple" << endreq;
      return StatusCode::FAILURE;
    }

    m_tupleBranch = m_tuple->Branch("XtalDigiTuple",
                                    &m_dat,
                                    m_tupleDesc);
    if (!m_tupleBranch) {
      msglog << MSG::ERROR << "Couldn't create tuple branch" << endreq;
      return StatusCode::FAILURE;
    }
  }

  return sc;
}

StatusCode XtalDigiTool::calculate(const CalXtalId &xtalId, 
                                   const vector<const Event::McIntegratingHit*> &hitList,
                                   const Event::EventHeader &evtHdr,            
                                   Event::CalDigi &calDigi,     
                                   bool &lacP,                  
                                   bool &lacN,                  
                                   bool &fleP,                  
                                   bool &fleN,                  
                                   bool &fheP,                  
                                   bool &fheN                   
                                   ) {
  StatusCode sc;

  m_dat.Clear();

  m_dat.twr = xtalId.getTower();
  m_dat.lyr = xtalId.getLayer();
  m_dat.col = xtalId.getColumn();

  /////////////////////////////////////
  // STAGE 0: retrieve calibrations ///
  /////////////////////////////////////
  sc = retrieveCalib(xtalId);
  if (sc.isFailure()) return sc;


  /////////////////////////////////////
  // STAGE 1: fill signal energies ////
  /////////////////////////////////////


  // loop over hits.
  for (vector<const Event::McIntegratingHit*>::const_iterator it = hitList.begin();
       it != hitList.end(); it++) {
    const Event::McIntegratingHit &hit = *(*it); // use ref to avoid ugly '->' operator

    // get volume identifier.
    VolumeIdentifier volId = ((VolumeIdentifier)hit.volumeID());

    // make sure the hits are cal hits
    if (!volId.isCal())
      throw invalid_argument("volume id is not Cal.  Programmer error");

    // make sure volumeid matches CalXtalId
    CalXtalId volXtalId(volId);
    if (volXtalId != xtalId)
      throw invalid_argument("volume id does not match xtalId.  Programmer error.");
    
    m_dat.sumEne += hit.totalEnergy();;

    ///////////////////////////////////////////////////////////////////////////
    ///////////////////// SUM DAC VALS FOR EACH HIT /////////////////////////
    ///////////////////////////////////////////////////////////////////////////

    //--XTAL DEPOSIT--//
    if((int)volId[fCellCmp] ==  m_eXtal ) {
      sc = sumCsIHit(xtalId, hit);
      if (sc.isFailure()) return sc;
    } 
    //--DIODE DEPOSIT--//
    else {
      sc = sumDiodeHit(xtalId, hit);
      if (sc.isFailure()) return sc;
    }
    m_dat.nMCHits++;
  }

  //////////////////////////////////////////////
  // STAGE 2:  add poissonic noise to signals //
  //////////////////////////////////////////////

  // -- poissonic noise is only added if we have hits.
  if (hitList.size() > 0 && !m_noRandNoise)
    for (XtalDiode xDiode; xDiode.isValid(); xDiode++) {
      DiodeNum diode = xDiode.getDiode();
      FaceNum  face  = xDiode.getFace();

      // convert dac in diode to MeV in xtal
      float meVXtal = m_dat.diodeDAC[face][diode]*m_dat.mpd[diode];

      // MeV in xtal -> electrons in diode
      float eDiode = meVXtal * m_ePerMeV[diode];

      // apply poissonic fluctuation to # of electrons.
      float noise = sqrt(eDiode)*RandGauss::shoot();

      // add noise
      eDiode += noise;

      // convert back to dac in diode
      meVXtal = eDiode/m_ePerMeV[diode];
      m_dat.diodeDAC[face][diode] = meVXtal/m_dat.mpd[diode];
    }     // for xDiode

  ///////////////////////////////
  // Stage 3: convert dac->adc //
  ///////////////////////////////

  for (FaceNum face; face.isValid(); face++)
    for (RngNum rng; rng.isValid(); rng++) {
      CalXtalId rngXtalId(m_dat.twr,m_dat.lyr,m_dat.col,
                          face, rng);
    
      // calcuate adc val from dac
      double tmp;
      sc = m_calCalibSvc->evalADC(rngXtalId, 
                                  m_dat.diodeDAC[face][rng.getDiode()],tmp);
      m_dat.adcPed[face][rng] = tmp;
      if (sc.isFailure()) return sc;   
    } // for xRng


  ///////////////////////////////
  // STAGE 4: electronic noise //
  ///////////////////////////////


  // electronic noise is calculated on a per-diode basis
  // uses sigmas from ADC ped data
  // adc vals need not include peds for this calculation
  if (!m_noRandNoise)
    for (XtalDiode xDiode; xDiode.isValid(); xDiode++) {
      DiodeNum diode = xDiode.getDiode();

      // only add electronic nose to sm diodes if we have
      // hits, otherwise may cause false signals?
      if (diode == SM_DIODE && !hitList.size()) continue;

      FaceNum face = xDiode.getFace();
    
      float sigX8 = m_dat.pedSig[face][diode.getX8Rng()];
      float sigX1 = m_dat.pedSig[face][diode.getX1Rng()];
    
      // use same rand for both since the X1 & X8 noise
      // is strongly correlated
      float rnd = RandGauss::shoot();

      m_dat.adcPed[face][diode.getX1Rng()] += sigX1*rnd;
      m_dat.adcPed[face][diode.getX8Rng()] += sigX8*rnd;    
    }


  ////////////////////////
  // Stage 5: LAC TESTS //
  ////////////////////////

  // test LACs against ped-subtracted adc values.
  if (m_allowNoiseOnlyTrig) { // noise CAN set lac trigger even if not hits.
    // set log-accept flags
    m_dat.lac[NEG_FACE] = m_dat.adcPed[NEG_FACE][LEX8] >= m_dat.lacThresh[NEG_FACE];
    m_dat.lac[POS_FACE] = m_dat.adcPed[POS_FACE][LEX8] >= m_dat.lacThresh[POS_FACE];
  }
  else 
    if (hitList.size() > 0) {
      m_dat.lac[NEG_FACE] = m_dat.adcPed[NEG_FACE][LEX8] >= m_dat.lacThresh[NEG_FACE];
      m_dat.lac[POS_FACE] = m_dat.adcPed[POS_FACE][LEX8] >= m_dat.lacThresh[POS_FACE];
    }


  ///////////////////////
  // Stage 6: Triggers //
  ///////////////////////
  sc = simTriggers();
  if (sc.isFailure()) return sc;
 
  
  ////////////////////////////////////////
  // Stage 7: ADD PEDS, CHECK ADC RANGE //
  ////////////////////////////////////////
  // double check that ADC vals are all >= 0
  // double check that ADC vals are all < maxadc(4095)
  for (XtalRng xRng; xRng.isValid(); xRng++) {
    FaceNum face = xRng.getFace();
    RngNum rng = xRng.getRng();

    m_dat.adc[face][rng] = m_dat.adcPed[face][rng] + m_dat.ped[face][rng];
    m_dat.adc[face][rng] = max<float>(0,m_dat.adc[face][rng]);
    m_dat.adc[face][rng] = min<float>(m_maxAdc,m_dat.adc[face][rng]);

    // round to nearest integer.
    m_dat.adc[face][rng] = (int)floor(m_dat.adc[face][rng]+0.5);
  }  // per range, add peds, check adc in range


    
  //////////////////////////////
  // Stage 8: Range Selection //
  //////////////////////////////
  sc = rangeSelect();
  if (sc.isFailure()) return sc;

  
  /////////////////////////////////////////////
  //-- STEP 9: Populate XtalDigiTuple vars --//
  /////////////////////////////////////////////
  
  if (m_tuple) {
    m_dat.RunID                   = evtHdr.run();
    m_dat.EventID                   = evtHdr.event();
    m_dat.success                = true;

    // calculate ene weighted pos
    if (m_dat.sumEneCsI) {// no divide-by-zero
      m_dat.csiWeightedPos /= m_dat.sumEneCsI;
    }
    
    // fill
    m_tuple->Fill();
  }

  
  ///////////////////////////////////////
  //-- STEP 10: Populate Return Vars --//
  ///////////////////////////////////////
  
  // generate xtalDigReadouts.
  sc = fillDigi(xtalId, calDigi);
  if (sc.isFailure()) return sc;

  lacP = (m_dat.lac[POS_FACE] != 0);
  lacN = (m_dat.lac[NEG_FACE] != 0);
  fleP = (m_dat.fle[POS_FACE] != 0);
  fleN = (m_dat.fle[NEG_FACE] != 0);
  fheP = (m_dat.fhe[POS_FACE] != 0);
  fheN = (m_dat.fhe[NEG_FACE] != 0);

  
  return StatusCode::SUCCESS;
}


StatusCode XtalDigiTool::retrieveCalib(const CalXtalId &xtalId) {
  StatusCode sc;

  //-- RETRIEVE MEV PER DAC--// 
  CalibData::ValSig mpdLrg, mpdSm;
  sc = m_calCalibSvc->getMeVPerDac(xtalId, mpdLrg, mpdSm);
  if (sc.isFailure()) return sc;
  m_dat.mpd[LRG_DIODE] = mpdLrg.getVal();
  m_dat.mpd[SM_DIODE]  = mpdSm.getVal();

  for (FaceNum face; face.isValid(); face++) {
    for (RngNum rng; rng.isValid(); rng++) {
      CalXtalId rngXtalId(m_dat.twr,m_dat.lyr,m_dat.col,face,rng);

      float tmp;
      //-- RETRIEVE PEDS --//
      sc = m_calCalibSvc->getPed(rngXtalId,
                                 m_dat.ped[face][rng],
                                 m_dat.pedSig[face][rng],
                                 tmp);
      if (sc.isFailure()) return sc;
      
      //-- RETRIEVE ULD --//
      CalibData::ValSig uldThold;
      sc = m_calCalibSvc->getULDCI(rngXtalId,uldThold);
      if (sc.isFailure()) return sc;
      m_dat.uldThold[face][rng] = uldThold.getVal();
    }
  
    CalXtalId faceXtalId(m_dat.twr,m_dat.lyr,m_dat.col,face);
    CalibData::ValSig fle,fhe,lac;
    sc = m_calCalibSvc->getTholdCI(faceXtalId,fle,fhe,lac);
    if (sc.isFailure()) return sc;
    m_dat.fleThresh[face] = fle.getVal();
    m_dat.fheThresh[face] = fhe.getVal();
    m_dat.lacThresh[face] = lac.getVal();
  }

  return StatusCode::SUCCESS;
  
}


StatusCode XtalDigiTool::sumCsIHit(const CalXtalId &xtalId, const Event::McIntegratingHit &hit) {
  StatusCode sc;

  float ene = hit.totalEnergy();
  
  float realpos = hit2pos(hit);

  //--ASYMMETRY--//
  double tmp;
  float asymL, asymS;
      
  sc = m_calCalibSvc->evalAsymLrg(xtalId, realpos, tmp);
  asymL = tmp;
  if (sc.isFailure()) return sc;
  sc = m_calCalibSvc->evalAsymSm(xtalId, realpos, tmp);
  asymS = tmp;
  if (sc.isFailure()) return sc;

  float meanDacS = ene/m_dat.mpd[SM_DIODE];
  float meanDacL = ene/m_dat.mpd[LRG_DIODE];

  //  calc dacSmP, dacSmN, dacLrgP, dacLrgN - here are quick notes
  //   asym=log(p/n)
  //   mean=sqrt(p*N)
  //   exp(asym)*mean^2=p/n*p*n=p^2
  //   p=exp(asym/2)*mean
  //   n=exp(-asym/2)*mean

  float dacSP = exp(   asymS/2) * meanDacS;
  float dacSN = exp(-1*asymS/2) * meanDacS;
  float dacLP = exp(   asymL/2) * meanDacL;
  float dacLN = exp(-1*asymL/2) * meanDacL;

  // sum energy to each diode
  m_dat.diodeDAC[POS_FACE][LRG_DIODE] += dacLP;
  m_dat.diodeDAC[POS_FACE][SM_DIODE]  += dacSP;
  m_dat.diodeDAC[NEG_FACE][LRG_DIODE] += dacLN;
  m_dat.diodeDAC[NEG_FACE][SM_DIODE]  += dacSN;

  // update summary vals
  if (m_tuple) {
    m_dat.nCsIHits++;
    m_dat.sumEneCsI += ene;
  
    // running total for ene weighted pos
    m_dat.csiWeightedPos += ene*realpos;
  }
  
  return StatusCode::SUCCESS;
}

float XtalDigiTool::hit2pos(const Event::McIntegratingHit &hit) {
  // get volume identifier.
  VolumeIdentifier volId = ((VolumeIdentifier)hit.volumeID());

  //-- HIT POSITION--//
  HepPoint3D mom1 = hit.moment1();
  int segm = volId[fSegment]; // segment # (0-11)
  // let's define the position of the segment along the crystal
  float relpos = (segm+0.5)/m_nCsISeg; // units in xtal len 0.0-1.0
  // in local reference system x is always oriented along the crystal
  float dpos =  mom1.x();

  // add position w/in segment to segment position
  relpos += dpos/m_CsILength; // still ranges 0.0-1.0

  // in 'mm' from center of xtal (units used in asym splines)
  return (relpos-0.5)*m_CsILength;
}

StatusCode XtalDigiTool::sumDiodeHit(const CalXtalId &xtalId, const Event::McIntegratingHit &hit) {
  StatusCode sc;
        
  float ene = hit.totalEnergy();
  float eDiode = ene*m_ePerMeVInDiode; // convert MeV-in-diode to electrons

  // get volume identifier.
  VolumeIdentifier volId = ((VolumeIdentifier)hit.volumeID());

  DiodeNum diode;
  FaceNum face;
  if((int)volId[fCellCmp]      == m_eDiodePLarge)  face = POS_FACE, diode = LRG_DIODE;
  else if((int)volId[fCellCmp] == m_eDiodeMLarge)  face = NEG_FACE, diode = LRG_DIODE;
  else if((int)volId[fCellCmp] == m_eDiodePSm)     face = POS_FACE, diode = SM_DIODE;
  else if((int)volId[fCellCmp] == m_eDiodeMSm)     face = NEG_FACE, diode = SM_DIODE;
  else throw invalid_argument("VolId is not in xtal or diode.  Programmer error.");
      
  // convert e-in-diode to MeV-in-xtal-center
  float meVXtal = eDiode/m_ePerMeV[diode]; 

  // convert MeV deposition to Dac val.
  // use asymmetry to determine how much of mean dac to apply to diode in question
  // remember we are treating mevXtal as at the xtal center (pos=0.0)
  float meanDAC;
  float asym;
  double tmp;
  meanDAC = meVXtal/m_dat.mpd[diode];
  if (diode == LRG_DIODE)
    sc = m_calCalibSvc->evalAsymLrg(xtalId,0.0,tmp);
  else
    sc = m_calCalibSvc->evalAsymSm(xtalId,0.0,tmp);
  asym = tmp;
  if (sc.isFailure()) return sc;

  // meanDAC = sqrt(dacP*dacN), exp(asym)=P/N, P=N*exp(asym), N=P/exp(asym)
  // mean^2 = P*N, mean^2=P*P/exp(asym), mean^2*exp(asym)=P^2
  float dacP;
  dacP = sqrt(meanDAC*meanDAC*exp(asym));
  float dac = (face==POS_FACE) ? dacP : dacP/exp(asym);

  // sum dac val to running total
  m_dat.diodeDAC[face][diode] += dac;

  m_dat.nDiodeHits[face][diode]++;

  return StatusCode::SUCCESS;
}

StatusCode XtalDigiTool::simTriggers() {
  StatusCode sc;
  
  for (FaceNum face; face.isValid(); face++) {

    //-- FLE --//
    float fleThresh = m_dat.fleThresh[face];

    // thresholds are in X8 units, but may need conversion to X1 units if they
    // are past hardware range.
    RngNum fleRng = (fleThresh > m_dat.uldThold[face][LEX8]) ?
      LEX1 : LEX8;
    
    // convert to LEX1 range if needed
    if (fleRng == LEX1) {
      //-- FIND RATIO LEX1/LEX8 by converting MAX LEX8 ADC ->LEX1 --//

      double x8ADC = m_dat.uldThold[face][LEX8];
      double tmpDAC;
      // 1st convert to dac
      sc = m_calCalibSvc->evalDAC(CalXtalId(m_dat.twr, m_dat.lyr, m_dat.col, face, LEX8),
                                  x8ADC, tmpDAC);
      if (sc.isFailure()) return sc;
      
      // 2nd convert to LEX1 adc
      double x1ADC;
      sc = m_calCalibSvc->evalADC(CalXtalId(m_dat.twr, m_dat.lyr, m_dat.col, face, LEX1),
                                  tmpDAC, x1ADC);
      if (sc.isFailure()) return sc;
      
      
      float rat = x1ADC/x8ADC;
      fleThresh *= rat;
    }

    // set trigger bit
    m_dat.fle[face] = (m_dat.adcPed[face][fleRng] >= fleThresh );


    //-- FHE --//
    float fheThresh = m_dat.fheThresh[face];

    // thresholds are in X8 units, but may need conversion to X1 units if they
    // are past hardware range.
    RngNum fheRng = (fheThresh > m_dat.uldThold[face][HEX8]) ?
      HEX1 : HEX8;
    
    // convert to HEX1 range if needed
    if (fheRng == HEX1) {
      //-- FIND RATIO HEX1/HEX8 by converting MAX HEX8 ADC -> HEX1 --//

      double x8ADC = m_dat.uldThold[face][HEX8];
      double tmpDAC;
      // 1st convert to dac
      sc = m_calCalibSvc->evalDAC(CalXtalId(m_dat.twr, m_dat.lyr, m_dat.col, face, HEX8),
                                  x8ADC, tmpDAC);
      if (sc.isFailure()) return sc;
      
      // 2nd convert to HEX1 adc
      double x1ADC;
      sc = m_calCalibSvc->evalADC(CalXtalId(m_dat.twr, m_dat.lyr, m_dat.col, face, HEX1),
                                  tmpDAC, x1ADC);
      if (sc.isFailure()) return sc;

      float rat = x1ADC/x8ADC;
      fheThresh *= rat;

    }

    // set trigger bit
    m_dat.fhe[face] = (m_dat.adcPed[face][fheRng] >= fheThresh);

    
  }

  return StatusCode::SUCCESS;
}


StatusCode XtalDigiTool::rangeSelect() {
  for (FaceNum face; face.isValid(); face++) {
    // BEST RANGE
    RngNum rng;
    for (rng=0; rng.isValid(); rng++) {
      // get ULD threshold
      CalXtalId tmpId(m_dat.twr, m_dat.lyr, m_dat.col,
                      face, rng);

      // case of HEX1 range 
      // adc vals have a ceiling in HEX1
      if (rng == HEX1) {
        if (m_dat.adcPed[face][rng] > m_dat.uldThold[face][HEX1]) {
          // set ADC to max val
          m_dat.adcPed[face][rng] =  m_dat.uldThold[face][HEX1];
          
          // set 'pegged' flag
          m_dat.saturated[face] = true;
        }
        
        // break before rng is incremented out-of-bounds
        break; 
      } else { // 1st 3 ranges
        // break on 1st energy rng that is < threshold
        if (m_dat.adcPed[face][rng] <= m_dat.uldThold[face][rng]) break;
      }
    }
    
    // assign range selection
    m_dat.rng[face] = (CalXtalId::AdcRange)(int)rng;
  }  // per face, range selection
  return StatusCode::SUCCESS;

}

StatusCode XtalDigiTool::fillDigi(const CalXtalId &xtalId, Event::CalDigi &calDigi) {
  //-- How many readouts ? --//
  int roLimit;
  short rangeMode = calDigi.getMode();
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

  // set status to ok for POS and NEG if no other bits set.
  unsigned short status = 0;
  // check for failure mode. If killed, set to zero and set DEAD bit
  if (m_FailSvc != 0) {  
    if (m_FailSvc->matchChannel(xtalId, (CalXtalId::POS)))
      if (m_dat.lac[POS_FACE]) (status = status | Event::CalDigi::CalXtalReadout::DEAD_P);
    if (m_FailSvc->matchChannel(xtalId, (CalXtalId::NEG)))
      if (m_dat.lac[NEG_FACE]) (status = status | Event::CalDigi::CalXtalReadout::DEAD_N);
  }

  if ((status & 0x00FF) == 0) status = 
                                (status | Event::CalDigi::CalXtalReadout::OK_P);
  if ((status & 0xFF00) == 0) status = 
                                (status | Event::CalDigi::CalXtalReadout::OK_N);

  // set up the digi
  for (int nRo=0; nRo < roLimit; nRo++) {
    // represents ranges used for current readout in loop
    short roRangeP = (m_dat.rng[POS_FACE] + nRo)%CalDefs::RngNum::N_VALS; 
    short roRangeN = (m_dat.rng[NEG_FACE] + nRo)%CalDefs::RngNum::N_VALS; 
          
    Event::CalDigi::CalXtalReadout ro = Event::CalDigi::CalXtalReadout(roRangeP, 
                                                                       (short)m_dat.adc[POS_FACE][roRangeP], 
                                                                       roRangeN, 
                                                                       (short)m_dat.adc[NEG_FACE][roRangeN], 
                                                                       status);
    calDigi.addReadout(ro);
  }


  return StatusCode::SUCCESS;
}

StatusCode XtalDigiTool::finalize() {
  // make sure optional tuple is closed out
  if (m_tupleFile) {
    m_tupleFile->Write();
    m_tupleFile->Close(); // trees deleted
    delete m_tupleFile;
    m_tupleFile = 0;
  }

  return StatusCode::SUCCESS;
}

