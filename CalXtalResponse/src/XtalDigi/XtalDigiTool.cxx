//    $Header$

/** @file
    @author Zach Fewtrell
*/

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

using namespace CalUtil;
using namespace Event;
using namespace CalXtalResponse;
using namespace CalibData;


static ToolFactory<XtalDigiTool> s_factory;
const IToolFactory& XtalDigiToolFactory = s_factory;

static float round_int(float in) { return floor(in + 0.5);}

XtalDigiTool::XtalDigiTool( const string& type,
                            const string& name,
                            const IInterface* parent)
  : AlgTool(type,name,parent),
    m_calCalibSvc(0),
    m_CsILength(0),
    m_eDiodeMSm(0),
    m_eDiodePSm(0),
    m_eDiodeMLarge(0),
    m_eDiodePLarge(0),
    m_ePerMeVInDiode(0),
    m_maxAdc(-1),
    m_tuple(0),
    m_tupleFile(0)
{
  declareInterface<IXtalDigiTool>(this);

  declareProperty("CalCalibSvc",        m_calCalibSvcName    = "CalCalibSvc");
  declareProperty("NoRandomNoise",      m_noRandNoise        = false);
  declareProperty("CalTrigTool",        m_calTrigToolName    = "CalTrigTool");
  declareProperty("tupleFilename",      m_tupleFilename      = "");
  declareProperty("AllowNoiseOnlyTrig", m_allowNoiseOnlyTrig = true);
  declareProperty("tupleLACOnly",       m_tupleLACOnly       = true);
  declareProperty("PrecalcCalibTool",   m_precalcCalibName   = "PrecalcCalibTool");
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
      msglog << MSG::ERROR << " constant " <<(*it).second <<" not defiPed" << endreq;
      return StatusCode::FAILURE;
    } else *((*it).first)=(int)val;

  // doubles are done separately
  double tmp;
  sc = detSvc->getNumericConstByName("CsILength", &tmp);
  m_CsILength = tmp;
  if (sc.isFailure()) {
    msglog << MSG::ERROR << " constant CsILength not defiPed" << endreq;
    return StatusCode::FAILURE;
  }
  
  // eperMeVInDiode originally from CalDigi XML file
  // but i don't want dependency, so i'm hard coding it.  (just a testing tool anyway :)
  m_ePerMeVInDiode = 2.77e5;
  
  // obtain CalCalibSvc
  sc = service(m_calCalibSvcName.value(), m_calCalibSvc);
  if (sc.isFailure()) {
    msglog << MSG::ERROR << "can't get " << m_calCalibSvcName  << endreq;
    return sc;
  }

  // open optional tuple file
  if (m_tupleFilename.value().length() > 0 ) {
    m_tupleFile = new TFile(m_tupleFilename.value().c_str(),"RECREATE","XtalDigiTuple");
    if (!m_tupleFile)
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
          !m_tuple->Branch("adcPed",         m_dat.adcPed.begin(),     "adcPed[2][4]/F")     ||
          !m_tuple->Branch("nMCHits",        &m_dat.nMCHits,           "nMCHits/i")          ||
          !m_tuple->Branch("nCsIHits",       &m_dat.nCsIHits,          "nCsIHits/i")         ||
          !m_tuple->Branch("nDiodeHits",     m_dat.nDiodeHits.begin(), "nDiodeHits[2][2]/i") ||
          !m_tuple->Branch("sumEneCsI",      &m_dat.sumEneCsI,         "sumEneCsI/F")        ||
          !m_tuple->Branch("sumEne",         &m_dat.sumEne,            "sumEne/F")           ||
          !m_tuple->Branch("diodeCIDAC",      m_dat.diodeCIDAC.begin(),"diodeCIDAC[2][2]/F") ||
          !m_tuple->Branch("csiWeightedPos", &m_dat.csiWeightedPos,    "csiWeightedPos/F")   ||
          !m_tuple->Branch("ped",            m_dat.ped.begin(),        "ped[2][4]/F")        ||
          !m_tuple->Branch("pedSigCIDAC",    m_dat.pedSigCIDAC.begin(),"pedSigCIDAC[2][4]/F")     ||
          !m_tuple->Branch("lacThreshCIDAC", m_dat.lacThreshCIDAC.begin(),  "lacThreshCIDAC[2]/F")     ||
          !m_tuple->Branch("trigThreshCIDAC", m_dat.trigThreshCIDAC.begin(), "trigThreshCIDAC[2][2]/F") ||
          !m_tuple->Branch("uldTholdADC",    m_dat.uldTholdADC.begin(),   "uldTholdADC[2][4]/F")   ||
          !m_tuple->Branch("mpd",            m_dat.mpd.begin(),        "mpd[2]/F")           ||
          !m_tuple->Branch("rng",            m_dat.rng.begin(),        "rng[2]/b")           ||
          !m_tuple->Branch("lac",            m_dat.lac.begin(),        "lac[2]/b")           ||
          !m_tuple->Branch("trigBits",       m_dat.trigBits.begin(),   "trigBits[2][2]/b")   ||
          !m_tuple->Branch("saturated",      m_dat.saturated.begin(),  "saturated[2]/b")) {
        msglog << MSG::ERROR << "Couldn't create tuple branch" << endreq;
        return StatusCode::FAILURE;
      }
    }
  } // optional tuple

  ///////////////////////////////////////
  //-- RETRIEVE HELPER TOOLS & SVCS  --//
  ///////////////////////////////////////

  // this tool needs to be shared by CalDigiAlg, XtalDigiTool & TriggerAlg, so I am
  // giving it global ownership
  sc = toolSvc()->retrieveTool("CalTrigTool", m_calTrigToolName, m_calTrigTool);
  if (sc.isFailure() ) {
    msglog << MSG::ERROR << "  Unable to create " << m_calTrigToolName << endreq;
    return sc;
  }
  
  // this tool may also be shared by CalTrigTool, global ownership
  sc = toolSvc()->retrieveTool("PrecalcCalibTool", m_precalcCalibName, m_precalcCalib);
  if (sc.isFailure() ) {
    msglog << MSG::ERROR << "  Unable to create " << m_precalcCalibName << endreq;
    return sc;
  }

  return StatusCode::SUCCESS;
}

/**
   - Loop through hitList & sum energies in CIDAC scale for each diode.
   - Apply direct diode deposits when specified by hit position.
   - Add noise, calculate adc output, range select, log accept & trigger output.
   - populate GltDigi class if glt != 0
*/
StatusCode XtalDigiTool::calculate(const vector<const McIntegratingHit*> &hitList,
                                   const EventHeader *evtHdr,            
                                   CalDigi &calDigi,     
                                   CalArray<FaceNum, bool> &lacBits,
                                   CalArray<XtalDiode, bool> &trigBits,
                                   Event::GltDigi *glt,
                                   bool zeroSuppress) {
  StatusCode sc;

  m_dat.Clear();

  CalXtalId xtalId = calDigi.getPackedId();
  m_dat.xtalIdx = XtalIdx(xtalId);

  /////////////////////////////////////
  // STAGE 1: fill signal energies ////
  /////////////////////////////////////
  
  // loop over hits.
  if (hitList.size()) {
    for (vector<const McIntegratingHit*>::const_iterator it = hitList.begin();
         it != hitList.end(); it++) {
      
      // only need MeVPerDAC if we have a deposit.
      const CalMevPerDac *mpd = m_calCalibSvc->getMPD(m_dat.xtalIdx);
      if (!mpd) return StatusCode::FAILURE;
      m_dat.mpd[LRG_DIODE] = mpd->getBig()->getVal();
      m_dat.mpd[SM_DIODE]  = mpd->getSmall()->getVal();
      
      const McIntegratingHit &hit = *(*it); // use ref to avoid ugly '->' operator
      
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
      ///////////////////// SUM CIDAC VALS FOR EACH HIT /////////////////////////
      ///////////////////////////////////////////////////////////////////////////
      
      //--XTAL DEPOSIT--//
      if((int)volId[fCellCmp] ==  m_eXtal ) {
        sc = sumCsIHit(hit);
        if (sc.isFailure()) return sc;
      } 
      //--DIODE DEPOSIT--//
      else {
        sc = sumDiodeHit(hit);
        if (sc.isFailure()) return sc;
      }
      m_dat.nMCHits++;
    } // per hit

    if (!m_noRandNoise)
      for (XtalDiode xDiode; xDiode.isValid(); xDiode++) {
        DiodeNum diode = xDiode.getDiode();
      
        //////////////////////////////////////////////
        // STAGE 2:  add poissonic noise to signals //
        //////////////////////////////////////////////
        // -- poissonic noise is only added if we have hits.
        
        // convert cidac in diode to MeV in xtal
        float meVXtal = m_dat.diodeCIDAC[xDiode]*m_dat.mpd[diode];
        
        // MeV in xtal -> electrons in diode
        float eDiode = meVXtal * m_ePerMeV[diode.val()];
        
        // apply poissonic fluctuation to # of electrons.
        float noise = sqrt(eDiode)*CLHEP::RandGauss::shoot();
        
        // add noise
        eDiode += noise;
        
        // convert back to cidac in diode
        meVXtal = eDiode/m_ePerMeV[diode.val()];
        m_dat.diodeCIDAC[xDiode] = 
          meVXtal/m_dat.mpd[diode];

      } // for (xDiode)
  } // if (hitList.size())  
    
  
  ///////////////////////////////
  // STAGE 3: electronic noise //
  ///////////////////////////////
  // electronic noise is calculated with sigmas from pedestal data
  // this calculation is required before LAC test & theoretically
  // should be done before trigger test as well.
  //
  // electronic noise is calculated in CIDAC units.  This allows
  // us to perform threhold checks before invoking dac2adc
  // spline functions.  In zero suppression mode, this should
  // allow a significant optimization.
  //
  // same noise is applied to both x1 && x8 range on single diode.
  for (FaceNum face; face.isValid(); face++) {
    if (!m_noRandNoise) 
      for (DiodeNum diode; diode.isValid(); diode++){
        // only add electronic nose to sm diodes if we have
        // hits, otherwise may cause false signals?
        if (diode == SM_DIODE && !hitList.size()) continue;
        
        XtalRng xRng(face, diode.getX8Rng());
        RngIdx rngIdx(m_dat.xtalIdx, xRng);
        
        // get pedestal sigma
        sc = m_precalcCalib->getPedSigCIDAC(rngIdx, m_dat.pedSigCIDAC[xRng]);
        if (sc.isFailure()) return sc;
        
        // use same rand for both ranges 
        // since the X1 & X8 noise
        float rnd = CLHEP::RandGauss::shoot();
        
        m_dat.diodeCIDAC[XtalDiode(face, diode)] += m_dat.pedSigCIDAC[xRng]*rnd;
      }
    
    ////////////////////////
    // Stage 4: LAC TESTS //
    ////////////////////////
    
    // test LACs against ped-subtracted adc values.
    if (m_allowNoiseOnlyTrig ||  // noise CAN set lac trigger even if not hits.
        hitList.size() > 0) {
      FaceIdx faceIdx(m_dat.xtalIdx, face);
      sc = m_precalcCalib->getLacCIDAC(faceIdx, m_dat.lacThreshCIDAC[face]);
      if (sc.isFailure()) return sc;
      
      // set log-accept flags
      lacBits[face] = m_dat.diodeCIDAC[XtalDiode(face, LRG_DIODE)] >= 
        m_dat.lacThreshCIDAC[face];
    }
  }

  ///////////////////////
  // Stage 5: Triggers //
  ///////////////////////
  sc = m_calTrigTool->calcXtalTrig(m_dat.xtalIdx,
                                   m_dat.diodeCIDAC,
                                   trigBits,
                                   glt);
  if (sc.isFailure()) return sc;

  //-- Quick Exit --//
  // once LAC & triggers have been processed (in CIDAC scale
  // I can now exit and save rest of processing including adc->dac
  // calculations
  if (zeroSuppress && !lacBits[POS_FACE] && !lacBits[NEG_FACE])
    return StatusCode::SUCCESS;

  
  ///////////////////////////////
  // Stage 6: convert cidac->adc //
  ///////////////////////////////
  for (XtalRng xRng; xRng.isValid(); xRng++) {
    XtalDiode xDiode(xRng.getFace(), xRng.getRng().getDiode());
    // calcuate adc val from cidac
    sc = m_calCalibSvc->evalADC(RngIdx(m_dat.xtalIdx, xRng), 
                                m_dat.diodeCIDAC[xDiode],
                                m_dat.adcPed[xRng]);
    if (sc.isFailure()) return sc;   
  } // xRng
  
  //////////////////////////////
  // Stage 7: Range Selection //
  //////////////////////////////
  sc = rangeSelect();
  if (sc.isFailure()) return sc;
  
  ////////////////////////////////////////////////////////
  //-- STEP 8: Populate XtalDigiTuple vars (optional) --//
  ////////////////////////////////////////////////////////
  
  // following steps only needed if tuple output is selected
  if (m_tuple) {
    // if lac_only, then only continue if one of 2 lac flags is high
    if ((m_tupleLACOnly && (lacBits[NEG_FACE] || 
                            lacBits[POS_FACE])) ||
        !m_tupleLACOnly) {
      
      if (evtHdr) {
        m_dat.RunID   = evtHdr->run();
        m_dat.EventID = evtHdr->event();
      }
      
      // calculate ene weighted pos
      if (m_dat.sumEneCsI)// no divide-by-zero
        m_dat.csiWeightedPos /= m_dat.sumEneCsI;
      
      copy(trigBits.begin(),
           trigBits.end(),
           m_dat.trigBits.begin());
      
      // instruct tuple to fill
      m_tuple->Fill();
    }
    
    // wack MSVC warning doesn't like me cast bool to int.
    m_dat.lac[POS_FACE] = lacBits[POS_FACE] != 0;
    m_dat.lac[NEG_FACE] = lacBits[NEG_FACE] != 0;
  }
  
  ///////////////////////////////////////
  //-- STEP 10: Populate Return Vars --//
  ///////////////////////////////////////
  // generate xtalDigReadouts.
  sc = fillDigi(calDigi);
  if (sc.isFailure()) return sc;
  
  return StatusCode::SUCCESS;
}

StatusCode XtalDigiTool::sumCsIHit(const Event::McIntegratingHit &hit) {
  StatusCode sc;

  float ene = hit.totalEnergy();
  
  float realpos = hit2pos(hit);

  //--ASYMMETRY--//
  float asymL, asymS;
      
  sc = m_calCalibSvc->evalAsym(m_dat.xtalIdx, ASYM_LL, realpos, asymL);
  if (sc.isFailure()) return sc;

  sc = m_calCalibSvc->evalAsym(m_dat.xtalIdx, ASYM_SS, realpos, asymS);
  if (sc.isFailure()) return sc;

  float meanCIDACS = ene/m_dat.mpd[SM_DIODE];
  float meanCIDACL = ene/m_dat.mpd[LRG_DIODE];

  //  calc cidacSmP, cidacSmN, cidacLrgP, cidacLrgN - here are quick notes
  //   asym=log(p/n)
  //   mean=sqrt(p*N)
  //   exp(asym)*mean^2=p/n*p*n=p^2
  //   p=exp(asym/2)*mean
  //   n=exp(-asym/2)*mean

  // sum energy to each diode
  m_dat.diodeCIDAC[XtalDiode(POS_FACE, LRG_DIODE)] += exp(   asymL/2) * meanCIDACL;
  m_dat.diodeCIDAC[XtalDiode(NEG_FACE, LRG_DIODE)] += exp(-1*asymL/2) * meanCIDACL;
  m_dat.diodeCIDAC[XtalDiode(POS_FACE, SM_DIODE)]  += exp(   asymS/2) * meanCIDACS;
  m_dat.diodeCIDAC[XtalDiode(NEG_FACE, SM_DIODE)]  += exp(-1*asymS/2) * meanCIDACS;

  // update summary vals
  if (m_tuple) {
    m_dat.nCsIHits++;
    m_dat.sumEneCsI += ene;
  
    // running total for ene weighted pos
    m_dat.csiWeightedPos += ene*realpos;
  }
  
  return StatusCode::SUCCESS;
}

float XtalDigiTool::hit2pos(const McIntegratingHit &hit) {
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

StatusCode XtalDigiTool::sumDiodeHit(const McIntegratingHit &hit) {
  StatusCode sc;
        
  float ene = hit.totalEnergy();
  float eDiode = ene*m_ePerMeVInDiode; // convert MeV-in-diode to electrons

  // get volume identifier.
  VolumeIdentifier volId = ((VolumeIdentifier)hit.volumeID());

  DiodeNum diode;
  FaceNum face;
  if((int)volId[fCellCmp]      == m_eDiodePLarge) face = POS_FACE, diode = LRG_DIODE;
  else if((int)volId[fCellCmp] == m_eDiodeMLarge) face = NEG_FACE, diode = LRG_DIODE;
  else if((int)volId[fCellCmp] == m_eDiodePSm)    face = POS_FACE, diode = SM_DIODE;
  else if((int)volId[fCellCmp] == m_eDiodeMSm)    face = NEG_FACE, diode = SM_DIODE;
  else throw invalid_argument("VolId is not in xtal or diode.  Programmer error.");
      
  // convert e-in-diode to MeV-in-xtal-center
  float meVXtal = eDiode/m_ePerMeV[diode.val()]; 

  // convert MeV deposition to CIDAC val.
  // use asymmetry to determine how much of mean cidac to apply to diode in question
  // remember we are treating mevXtal as at the xtal center (pos=0.0)
  float meanCIDAC;
  float asym;
  meanCIDAC = meVXtal/m_dat.mpd[diode];
  if (diode == LRG_DIODE)
    sc = m_calCalibSvc->getAsymCtr(m_dat.xtalIdx, ASYM_LL, asym);
  else
    sc = m_calCalibSvc->getAsymCtr(m_dat.xtalIdx, ASYM_SS, asym);
  if (sc.isFailure()) return sc;

  // meanCIDAC = sqrt(cidacP*cidacN), exp(asym)=P/N, P=N*exp(asym), N=P/exp(asym)
  // mean^2 = P*N, mean^2=P*P/exp(asym), mean^2*exp(asym)=P^2
  float cidacP;
  cidacP = sqrt(meanCIDAC*meanCIDAC*exp(asym));
  float cidac = (face==POS_FACE) ? cidacP : cidacP/exp(asym);

  // sum cidac val to running total
  XtalDiode xDiode(face,diode);
  m_dat.diodeCIDAC[xDiode] += cidac;
  m_dat.nDiodeHits[xDiode]++;

  return StatusCode::SUCCESS;
}

StatusCode XtalDigiTool::rangeSelect() {
  for (FaceNum face; face.isValid(); face++) {
    // retrieve calibration

    //-- TholdCI --//
    FaceIdx faceIdx(m_dat.xtalIdx, face);
    const CalTholdCI *tholdCI = m_calCalibSvc->getTholdCI(faceIdx);
    if (!tholdCI) return StatusCode::FAILURE;;

    // BEST RANGE
    RngNum rng;
    for (rng=0; rng.isValid(); rng++) {
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
    m_dat.rng[face] = (CalXtalId::AdcRange)rng;
  }  // per face, range selection

  return StatusCode::SUCCESS;
}

StatusCode XtalDigiTool::fillDigi(CalDigi &calDigi) {
  
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

  // set up the digi
  CalArray<FaceNum, RngNum> roRange;
  CalArray<FaceNum, float> adc;
  for (int nRo=0; nRo < roLimit; nRo++) {
    for (FaceNum face; face.isValid(); face++) {
      // represents ranges used for current readout in loop
      roRange[face] = (m_dat.rng[face].val() + nRo) % RngNum::N_VALS; 

      XtalRng xRng(face,roRange[face]);
	  RngIdx rngIdx(m_dat.xtalIdx, face, roRange[face]);
	  const CalibData::Ped *ped = m_calCalibSvc->getPed(rngIdx);
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
                                                         0);
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

  
  
