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


static ToolFactory<XtalDigiTool> s_factory;
const IToolFactory& XtalDigiToolFactory = s_factory;

static float rint(float in) { return floor(in + 0.5);}

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
    m_tupleFile(0)
{
  declareInterface<IXtalDigiTool>(this);

  declareProperty("CalCalibSvc",        m_calCalibSvcName    = "CalCalibSvc");
  declareProperty("NoRandomNoise",      m_noRandNoise        = false);
  declareProperty("CalTrigTool",        m_calTrigToolName    = "CalTrigTool");
  declareProperty("tupleFilename",      m_tupleFilename      = "");
  declareProperty("AllowNoiseOnlyTrig", m_allowNoiseOnlyTrig = true);
  declareProperty("tupleLACOnly",       m_tupleLACOnly       = true);
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
          !m_tuple->Branch("adc",            m_dat.adc.begin(),        "adc[2][4]/F")        ||
          !m_tuple->Branch("adcPed",         m_dat.adcPed.begin(),     "adcPed[2][4]/F")     ||
          !m_tuple->Branch("nMCHits",        &m_dat.nMCHits,           "nMCHits/i")          ||
          !m_tuple->Branch("nCsIHits",       &m_dat.nCsIHits,          "nCsIHits/i")         ||
          !m_tuple->Branch("nDiodeHits",     m_dat.nDiodeHits.begin(), "nDiodeHits[2][2]/i") ||
          !m_tuple->Branch("sumEneCsI",      &m_dat.sumEneCsI,         "sumEneCsI/F")        ||
          !m_tuple->Branch("sumEne",         &m_dat.sumEne,            "sumEne/F")           ||
          !m_tuple->Branch("diodeDAC",       m_dat.diodeDAC.begin(),   "diodeDAC[2][2]/F")   ||
          !m_tuple->Branch("csiWeightedPos", &m_dat.csiWeightedPos,    "csiWeightedPos/F")   ||
          !m_tuple->Branch("ped",            m_dat.ped.begin(),        "ped[2][4]/F")        ||
          !m_tuple->Branch("pedSig",         m_dat.pedSig.begin(),     "pedSig[2][4]/F")     ||
          !m_tuple->Branch("lacThresh",      m_dat.lacThresh.begin(),  "lacThresh[2]/F")     ||
          !m_tuple->Branch("trigThresh",     m_dat.trigThresh.begin(), "trigThresh[2][2]/F") ||
          !m_tuple->Branch("uldThold",       m_dat.uldThold.begin(),   "uldThold[2][4]/F")   ||
          !m_tuple->Branch("mpd",            m_dat.mpd.begin(),        "mpd[2]/F")           ||
          !m_tuple->Branch("twr",            &m_dat.twr,               "twr/b")              ||
          !m_tuple->Branch("lyr",            &m_dat.lyr,               "lyr/b")              ||
          !m_tuple->Branch("col",            &m_dat.col,               "col/b")              ||
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

  return StatusCode::SUCCESS;
}

/**
   - Loop through hitList & sum energies in DAC scale for each diode.
   - Apply direct diode deposits when specified by hit position.
   - Add noise, calculate adc output, range select, log accept & trigger output.
   - populate GltDigi class if glt != 0
*/
StatusCode XtalDigiTool::calculate(const vector<const McIntegratingHit*> &hitList,
                                   const EventHeader *evtHdr,            
                                   CalDigi &calDigi,     
                                   CalArray<FaceNum, bool> &lacBits,
                                   CalArray<XtalDiode, bool> &trigBits,
                                   Event::GltDigi *glt) {
  StatusCode sc;

  m_dat.Clear();

  m_dat.xtalIdx = XtalIdx(calDigi.getPackedId());
  m_dat.twr = m_dat.xtalIdx.getTwr();
  m_dat.lyr = m_dat.xtalIdx.getLyr();
  m_dat.col = m_dat.xtalIdx.getCol();

  /////////////////////////////////////
  // STAGE 0: retrieve calibrations ///
  /////////////////////////////////////
  sc = retrieveCalib();
  if (sc.isFailure()) return sc;

  // only need MeVPerDAC if we have a deposit.
  if (hitList.size()) 
    sc = m_calCalibSvc->getMPD(m_dat.xtalIdx, m_dat.mpd);
  if (sc.isFailure()) return sc;


  /////////////////////////////////////
  // STAGE 1: fill signal energies ////
  /////////////////////////////////////


  // loop over hits.
  for (vector<const McIntegratingHit*>::const_iterator it = hitList.begin();
       it != hitList.end(); it++) {
    const McIntegratingHit &hit = *(*it); // use ref to avoid ugly '->' operator

    // get volume identifier.
    VolumeIdentifier volId = ((VolumeIdentifier)hit.volumeID());

    // make sure the hits are cal hits
    if (!volId.isCal())
      throw invalid_argument("volume id is not Cal.  Programmer error");

    // make sure volumeid matches CalXtalId
    CalXtalId volXtalId(volId);
    if (volXtalId != m_dat.xtalIdx.getCalXtalId())
      throw invalid_argument("volume id does not match xtalId.  Programmer error.");
    
    m_dat.sumEne += hit.totalEnergy();;

    ///////////////////////////////////////////////////////////////////////////
    ///////////////////// SUM DAC VALS FOR EACH HIT /////////////////////////
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

  for (FaceNum face; face.isValid(); face++) {
    for (DiodeNum diode; diode.isValid(); diode++) {
      XtalDiode xDiode(face,diode);

      //////////////////////////////////////////////
      // STAGE 2:  add poissonic noise to signals //
      //////////////////////////////////////////////
      // -- poissonic noise is only added if we have hits.
      if (hitList.size() > 0 && !m_noRandNoise) {
        // convert dac in diode to MeV in xtal
        float meVXtal = m_dat.diodeDAC[xDiode]*m_dat.mpd[diode];
        
        // MeV in xtal -> electrons in diode
        float eDiode = meVXtal * m_ePerMeV[diode.getInt()];
        
        // apply poissonic fluctuation to # of electrons.
        float noise = sqrt(eDiode)*RandGauss::shoot();
        
        // add noise
        eDiode += noise;
        
        // convert back to dac in diode
        meVXtal = eDiode/m_ePerMeV[diode.getInt()];
        m_dat.diodeDAC[xDiode] = meVXtal/m_dat.mpd[diode];
      } // poissonic noise
        

      ///////////////////////////////
      // Stage 3: convert dac->adc //
      ///////////////////////////////
      for (THXNum thx; thx.isValid(); thx++) {
        RngNum rng(diode,thx);
        RngIdx rngIdx(m_dat.twr,m_dat.lyr,m_dat.col,
                      face, rng);
                  
        XtalRng xRng(face,rng);          
        // calcuate adc val from dac
        sc = m_calCalibSvc->evalADC(rngIdx, 
                                    m_dat.diodeDAC[xDiode],
                                    m_dat.adcPed[xRng]);
        if (sc.isFailure()) return sc;   

        ///////////////////////////////
        // STAGE 4: electronic noise //
        ///////////////////////////////
        // electronic noise is calculated on a per-diode basis
        // uses sigmas from ADC ped data
        // adc vals need not include peds for this calculation
        if (!m_noRandNoise) {
          // only add electronic nose to sm diodes if we have
          // hits, otherwise may cause false signals?
          if (diode == SM_DIODE && !hitList.size()) continue;
        
          float sig = m_dat.pedSig[xRng];
          
          // use same rand for both since the X1 & X8 noise
          // is strongly correlated
          float rnd = RandGauss::shoot();
        
          m_dat.adcPed[xRng] += sig*rnd;
        }


      }  // thx (x8,x1)

    } // diode
    
    ////////////////////////
    // Stage 6: LAC TESTS //
    ////////////////////////
    
    // test LACs against ped-subtracted adc values.
    if (m_allowNoiseOnlyTrig ||  // noise CAN set lac trigger even if not hits.
        hitList.size() > 0) {
      // set log-accept flags
      m_dat.lac[face] = m_dat.adcPed[XtalRng(face, LEX8)] >= 
        m_dat.lacThresh[face];
    }

  } // face loop

  ///////////////////////
  // Stage 7: Triggers //
  ///////////////////////
  sc = m_calTrigTool->calcXtalTrig(m_dat.xtalIdx,
                                   m_dat.adcPed,
                                   trigBits,
                                   glt);
  if (sc.isFailure()) return sc;

    
  //////////////////////////////
  // Stage 8: Range Selection //
  //////////////////////////////
  sc = rangeSelect();
  if (sc.isFailure()) return sc;

  
  ////////////////////////////////////////
  // Stage 5: ADD PEDS, CHECK ADC RANGE //
  ////////////////////////////////////////
  // double check that ADC vals are all >= 0
  // double check that ADC vals are all < maxadc(4095)

  // must be after rangeSelect()  b/c rangeSelect()
  // may clip HEX1 to HEX1 saturation point
  for (XtalRng xRng; xRng.isValid(); xRng++) {
    m_dat.adc[xRng] = 
      m_dat.adcPed[xRng] + m_dat.ped[xRng];
    m_dat.adc[xRng] = 
      max<float>(0,m_dat.adc[xRng]);
    m_dat.adc[xRng] = 
      min<float>(m_maxAdc,m_dat.adc[xRng]);

    //-- round adc values to nearest int.
    m_dat.adc[xRng] = rint(m_dat.adc[xRng]);
  }

  
  /////////////////////////////////////////////
  //-- STEP 9: Populate XtalDigiTuple vars --//
  /////////////////////////////////////////////
  
  // following steps only needed if tuple output is selected
  if (m_tuple) {

    // if lac_only, then only continue if one of 2 lac flags is high
    if ((m_tupleLACOnly && (m_dat.lac[NEG_FACE] || 
                            m_dat.lac[POS_FACE])) ||
        !m_tupleLACOnly) {
         
      
      if (evtHdr) {
        m_dat.RunID   = evtHdr->run();
        m_dat.EventID = evtHdr->event();
      }

      // calculate ene weighted pos
      if (m_dat.sumEneCsI) {// no divide-by-zero
        m_dat.csiWeightedPos /= m_dat.sumEneCsI;
      }

      copy(trigBits.begin(),
           trigBits.end(),
           m_dat.trigBits.begin());
    
      // instruct tuple to fill
      m_tuple->Fill();
    }
  }

  
  ///////////////////////////////////////
  //-- STEP 10: Populate Return Vars --//
  ///////////////////////////////////////


  // generate xtalDigReadouts.
  sc = fillDigi(calDigi);
  if (sc.isFailure()) return sc;

  copy(m_dat.lac.begin(), 
       m_dat.lac.end(), 
       lacBits.begin());

  return StatusCode::SUCCESS;
}

StatusCode XtalDigiTool::retrieveCalib() {
  StatusCode sc;

  //-- RETRIEVE PEDS --//
  sc = m_calCalibSvc->getPed(m_dat.xtalIdx,
                             m_dat.ped,
                             m_dat.pedSig);
  if (sc.isFailure()) return sc;

  //-- RETRIEVE ULD --//
  sc = m_calCalibSvc->getULDCI(m_dat.xtalIdx,
                               m_dat.uldThold);
  if (sc.isFailure()) return sc;
                               

  //-- RETRIEVE THOLDCI --//
  sc = m_calCalibSvc->getTholdCI(m_dat.xtalIdx,
                                 m_dat.trigThresh,
                                 m_dat.lacThresh);
  if (sc.isFailure()) return sc;
  
  return StatusCode::SUCCESS;
}


StatusCode XtalDigiTool::sumCsIHit(const McIntegratingHit &hit) {
  StatusCode sc;

  float ene = hit.totalEnergy();
  
  float realpos = hit2pos(hit);

  //--ASYMMETRY--//
  float asymL, asymS;
      
  sc = m_calCalibSvc->evalAsymLrg(m_dat.xtalIdx, realpos, asymL);
  if (sc.isFailure()) return sc;

  sc = m_calCalibSvc->evalAsymSm(m_dat.xtalIdx, realpos, asymS);
  if (sc.isFailure()) return sc;

  float meanDacS = ene/m_dat.mpd[SM_DIODE];
  float meanDacL = ene/m_dat.mpd[LRG_DIODE];

  //  calc dacSmP, dacSmN, dacLrgP, dacLrgN - here are quick notes
  //   asym=log(p/n)
  //   mean=sqrt(p*N)
  //   exp(asym)*mean^2=p/n*p*n=p^2
  //   p=exp(asym/2)*mean
  //   n=exp(-asym/2)*mean

  // sum energy to each diode
  m_dat.diodeDAC[XtalDiode(POS_FACE, LRG_DIODE)] += exp(   asymL/2) * meanDacL;
  m_dat.diodeDAC[XtalDiode(NEG_FACE, LRG_DIODE)] += exp(-1*asymL/2) * meanDacL;
  m_dat.diodeDAC[XtalDiode(POS_FACE, SM_DIODE)]  += exp(   asymS/2) * meanDacS;
  m_dat.diodeDAC[XtalDiode(NEG_FACE, SM_DIODE)]  += exp(-1*asymS/2) * meanDacS;

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
  float meVXtal = eDiode/m_ePerMeV[diode.getInt()]; 

  // convert MeV deposition to Dac val.
  // use asymmetry to determine how much of mean dac to apply to diode in question
  // remember we are treating mevXtal as at the xtal center (pos=0.0)
  float meanDAC;
  float asym;
  meanDAC = meVXtal/m_dat.mpd[diode];
  if (diode == LRG_DIODE)
    sc = m_calCalibSvc->evalAsymLrg(m_dat.xtalIdx,0.0,asym);
  else
    sc = m_calCalibSvc->evalAsymSm(m_dat.xtalIdx,0.0,asym);
  if (sc.isFailure()) return sc;

  // meanDAC = sqrt(dacP*dacN), exp(asym)=P/N, P=N*exp(asym), N=P/exp(asym)
  // mean^2 = P*N, mean^2=P*P/exp(asym), mean^2*exp(asym)=P^2
  float dacP;
  dacP = sqrt(meanDAC*meanDAC*exp(asym));
  float dac = (face==POS_FACE) ? dacP : dacP/exp(asym);

  // sum dac val to running total
  XtalDiode xDiode(face,diode);
  m_dat.diodeDAC[xDiode] += dac;
  m_dat.nDiodeHits[xDiode]++;

  return StatusCode::SUCCESS;
}

StatusCode XtalDigiTool::rangeSelect() {
  for (FaceNum face; face.isValid(); face++) {
    // BEST RANGE
    RngNum rng;
    for (rng=0; rng.isValid(); rng++) {
      // get ULD threshold
      XtalRng xRng(face,rng);

      // case of HEX1 range 
      // adc vals have a ceiling in HEX1
      if (rng == HEX1) {
        if (m_dat.adcPed[xRng] >= m_dat.uldThold[xRng]) {
          // set ADC to max val
          m_dat.adcPed[xRng] =  m_dat.uldThold[xRng];
          
          // set 'pegged' flag
          m_dat.saturated[face] = true;
        }
        
        // break before rng is incremented out-of-bounds
        break; 
      } else { // 1st 3 ranges
        // break on 1st energy rng that is < threshold
        if (m_dat.adcPed[xRng] <= m_dat.uldThold[xRng]) break;
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

  // set status to ok for POS and NEG if no other bits set.
  unsigned short status = 0;
  // check for failure mode. If killed, set to zero and set DEAD bit
  if (m_FailSvc != 0) {  
    if (m_FailSvc->matchChannel(m_dat.xtalIdx.getCalXtalId(), (CalXtalId::POS)))
      if (m_dat.lac[POS_FACE]) (status = status | CalDigi::CalXtalReadout::DEAD_P);
    if (m_FailSvc->matchChannel(m_dat.xtalIdx.getCalXtalId(), (CalXtalId::NEG)))
      if (m_dat.lac[NEG_FACE]) (status = status | CalDigi::CalXtalReadout::DEAD_N);
  }

  if ((status & 0x00FF) == 0) status = 
                                (status | CalDigi::CalXtalReadout::OK_P);
  if ((status & 0xFF00) == 0) status = 
                                (status | CalDigi::CalXtalReadout::OK_N);

  // set up the digi
  for (int nRo=0; nRo < roLimit; nRo++) {
    // represents ranges used for current readout in loop
    short roRangeP = (m_dat.rng[POS_FACE].getInt() + nRo)%RngNum::N_VALS; 
    short roRangeN = (m_dat.rng[NEG_FACE].getInt() + nRo)%RngNum::N_VALS; 
          
    CalDigi::CalXtalReadout ro = CalDigi::CalXtalReadout(roRangeP, 
                                                         (short)m_dat.adc[XtalRng(POS_FACE, roRangeP)], 
                                                         roRangeN, 
                                                         (short)m_dat.adc[XtalRng(NEG_FACE, roRangeN)], 
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

