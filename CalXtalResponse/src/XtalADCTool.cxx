// Include files

// LOCAL
#include "CalXtalResponse/IXtalADCTool.h"
#include "CalXtalResponse/ICalCalibSvc.h"

// GLAST
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "CalUtil/CalDefs.h"

// EXTLIB
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/MsgStream.h"

#include "CLHEP/Random/RandGauss.h"

// STD

using namespace CalDefs;

/*! \class XtalADCTool
  \author Zachary Fewtrell
  \brief Official implementation of IXtalADCTool.  sums IntMcHits into digital response for one Cal xtalId.

*/

class XtalADCTool : public AlgTool, virtual public IXtalADCTool {
public:
  /// default ctor, declares jobOptions
  XtalADCTool::XtalADCTool( const string& type, 
                            const string& name, 
                            const IInterface* parent);

  /// retrieves needed parameters and pointers to required services
  virtual StatusCode initialize();

  /// calculate xtal Adc response for all rngs basedon collection of McHits in xtal & diode regions
  StatusCode calculate(const CalXtalId &xtalId, 
                       const vector<const Event::McIntegratingHit*> &hitList, 
                       bool &lacP,
                       bool &lacN,
                       CalXtalId::AdcRange &rngP,
                       CalXtalId::AdcRange &rngN,
                       vector<int> &adcP,                
                       vector<int> &adcN,                
                       bool &peggedP,                     
                       bool &peggedN
                       );
private:
  StringProperty m_calCalibSvcName;      ///< name of CalCalibSvc to use for calib constants.
  ICalCalibSvc *m_calCalibSvc;           ///< pointer to CalCalibSvc object.

  //-- Gaudi supplId constants --//
  int m_nCsISeg;                         ///< number of geometric segments per Xtal
  int m_eXtal;
  double m_CsILength;                    ///< Xtal length
  int m_eDiodeMSm;                    ///< detModel identifier for sm minus-side diode
  int m_eDiodePSm;                    ///< detModel identifier for sm plus-side diode
  int m_eDiodeMLarge;                    ///< detModel identifier for large minus-side diode
  int m_eDiodePLarge;                    ///< detModel identifier for large plus-side diode
  double m_ePerMeVInDiode;
  int m_ePerMeV[2];                      ///< gain - electrons/MeV 1=Sm, 0=Large
  int m_maxAdc;                          ///< max value for ADC
};

static ToolFactory<XtalADCTool> s_factory;
const IToolFactory& XtalADCToolFactory = s_factory;

XtalADCTool::XtalADCTool( const string& type,
                          const string& name,
                          const IInterface* parent)
  : AlgTool(type,name,parent) {
  declareInterface<IXtalADCTool>(this);

  declareProperty("CalCalibSvc", m_calCalibSvcName="CalCalibSvc");
}

StatusCode XtalADCTool::initialize() {
  MsgStream msglog(msgSvc(), name());
  msglog << MSG::INFO << "initialize" << endreq;

  StatusCode sc;

  double value;
  typedef map<int*,string> PARAMAP;

  // try to find the GlastDevSvc service
  IGlastDetSvc* detSvc;
  sc = service("GlastDetSvc", detSvc);
  if (sc.isFailure() ) {
    msglog << MSG::ERROR << "  Unable to retrieve GlastDetSvc " << endreq;
    return sc;
  }

  //-- Retrieve GlastDetSvc constants --//
  PARAMAP param;
  param[&m_nCsISeg]     = string("nCsISeg");
  param[&m_eXtal]       = string("eXtal");
  param[&m_eDiodeMSm]= string("eDiodeMSmall");
  param[&m_eDiodePSm]= string("eDiodePSmall");
  param[&m_eDiodeMLarge]= string("eDiodeMLarge");
  param[&m_eDiodePLarge]= string("eDiodePLarge");
  param[m_ePerMeV+1]    = string("cal.ePerMeVSmall");
  param[m_ePerMeV]      = string("cal.ePerMevLarge");
  param[&m_maxAdc]      = string("cal.maxAdcValue");

  for(PARAMAP::iterator it=param.begin(); it!=param.end();it++)
    if(!detSvc->getNumericConstByName((*it).second, &value)) {
      msglog << MSG::ERROR << " constant " <<(*it).second <<" not defined" << endreq;
      return StatusCode::FAILURE;
    } else *((*it).first)=(int)value;

  typedef map<double*,string> DPARAMAP;
  DPARAMAP dparam;
  dparam[&m_CsILength]  = string("CsILength");

  for(DPARAMAP::iterator dit=dparam.begin(); dit!=dparam.end();dit++)
    if(!detSvc->getNumericConstByName((*dit).second,(*dit).first)) {
      msglog << MSG::ERROR << " constant " <<(*dit).second << " not defined" << endreq;
      return StatusCode::FAILURE;
    }
  
  // eperMeVInDiode originally from CalDigi XML file
  // but i don't want dependency, so i'm hard coding it.  (just a testing tool anyway :)
  m_ePerMeVInDiode = 2.77e5;
  
  // obtain CalCalibSvc
  sc = service(m_calCalibSvcName.value(), m_calCalibSvc);
  if (sc.isFailure()) {
    msglog << MSG::ERROR << "Unable to get CalCalibSvc." << endreq;
    return sc;
  }

  return sc;
}

StatusCode XtalADCTool::calculate(const CalXtalId &xtalId,
                                  const vector<const Event::McIntegratingHit*> &hitList, // list of all mc hits for this xtal & it's diodes.
                                  bool &lacP,
                                  bool &lacN,
                                  CalXtalId::AdcRange &rngP,
                                  CalXtalId::AdcRange &rngN,
                                  vector<int> &adcP,
                                  vector<int> &adcN,
                                  bool &peggedP,
                                  bool &peggedN
                                  ) {
  MsgStream msglog(msgSvc(), name());
  StatusCode sc;

  CalVec<XtalDiode, double> diodeDAC(XtalDiode::N_VALS);
  CalVec<XtalRng, double> tmpADC(XtalRng::N_VALS);

  // STAGE 1: fill signal energies /////////////////////////////////////////////////////

  //--MEV PER DAC--// - will be used throughout the algorithm
  CalibData::ValSig mpdLrg, mpdSm;
  sc = m_calCalibSvc->getMeVPerDac(xtalId, mpdLrg, mpdSm);
  if (sc.isFailure()) return sc;

  if (hitList.size())
     msglog << MSG::DEBUG << "xtalId=" << xtalId << "\tnHits=" << hitList.size() << endreq;
    
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

    double ene = hit.totalEnergy();
    msglog << MSG::DEBUG << "\tcell=" << volId[fCellCmp] << " ene=" << ene << endreq;

    ///////////////////////////////////////////////////////////////////////////
    ///////////////////// SUM DAC VALUES FOR EACH HIT /////////////////////////
    ///////////////////////////////////////////////////////////////////////////

    //--XTAL DEPOSIT--//
    
    if((int)volId[fCellCmp] ==  m_eXtal ) {
      //-- HIT POSITION--//
      HepPoint3D mom1 = hit.moment1();
      int segm = volId[fSegment]; // segment # (0-11)
      // let's define the position of the segment along the crystal
      double relpos = (segm+0.5)/m_nCsISeg; // units in xtal len 0.0-1.0
      // in local reference system x is always oriented along the crystal
      double dpos =  mom1.x();

      // add position w/in segment to segment position
      relpos += dpos/m_CsILength; // still ranges 0.0-1.0

      // in 'mm' from center of xtal (units used in asym splines)
      double realpos = (relpos-0.5)*m_CsILength;

      //--ASYMMETRY--//
      double asymL, asymS;
      
      sc = m_calCalibSvc->evalAsymLrg(xtalId, realpos, asymL);
      if (sc.isFailure()) return sc;
      sc = m_calCalibSvc->evalAsymSm(xtalId, realpos, asymS);
      if (sc.isFailure()) return sc;

      double meanDacS = ene/mpdSm.getVal();
      double meanDacL = ene/mpdLrg.getVal();

      //  calc dacSmP, dacSmN, dacLrgP, dacLrgN - here are quick notes
      //   asym=log(p/n)
      //   mean=sqrt(p*N)
      //   exp(asym)*mean^2=p/n*p*n=p^2
      //   p=exp(asym/2)*mean
      //   n=exp(-asym/2)*mean

      double dacSP = exp(   asymS/2) * meanDacS;
      double dacSN = exp(-1*asymS/2) * meanDacS;
      double dacLP = exp(   asymL/2) * meanDacL;
      double dacLN = exp(-1*asymL/2) * meanDacL;

      // sum energy to each diode
      diodeDAC[XtalDiode(POS_FACE,LRG_DIODE)] += dacLP;
      diodeDAC[XtalDiode(POS_FACE,SM_DIODE)]  += dacSP;
      diodeDAC[XtalDiode(NEG_FACE,LRG_DIODE)] += dacLN;
      diodeDAC[XtalDiode(NEG_FACE,SM_DIODE)]  += dacSN;
    } else {
    
      //--DIRECT DIODE DEPOSIT--//

      double eDiode = ene*m_ePerMeVInDiode; // convert MeV-in-diode to electrons

      DiodeNum diode;
      FaceNum face;
      if((int)volId[fCellCmp]      == m_eDiodePLarge)  face = POS_FACE, diode = LRG_DIODE;
      else if((int)volId[fCellCmp] == m_eDiodeMLarge)  face = NEG_FACE, diode = LRG_DIODE;
      else if((int)volId[fCellCmp] == m_eDiodePSm)     face = POS_FACE, diode = SM_DIODE;
      else if((int)volId[fCellCmp] == m_eDiodeMSm)     face = NEG_FACE, diode = SM_DIODE;
      else throw invalid_argument("VolId is not in xtal or diode.  Programmer error.");
      
      // convert e-in-diode to MeV-in-xtal-center
      double meVXtal = eDiode/m_ePerMeV[diode]; 

      // convert MeV deposition to Dac value.
      // use asymmetry to determine how much of mean dac to apply to diode in question
      // remember we are treating mevXtal as at the xtal center (pos=0.0)
      double meanDAC;
      double asym;
      if (diode == LRG_DIODE) {
        meanDAC = meVXtal/mpdLrg.getVal();
        sc = m_calCalibSvc->evalAsymLrg(xtalId,0.0,asym);
      } else {
        meanDAC = meVXtal/mpdSm.getVal();
        sc = m_calCalibSvc->evalAsymSm(xtalId,0.0,asym);
      }
      if (sc.isFailure()) return sc;

      // meanDAC = sqrt(dacP*dacN), exp(asym)=P/N, P=N*exp(asym), N=P/exp(asym)
      // mean^2 = P*N, mean^2=P*P/exp(asym), mean^2*exp(asym)=P^2
      double dacP;
      dacP = sqrt(meanDAC*meanDAC*exp(asym));
      double dac = (face=POS_FACE) ? dacP : dacP/exp(asym);

      // sum dac value to running total
      diodeDAC[XtalDiode(face,diode)] += dac;
    }
  }

  // STAGE 2:  add poissonic noise to signals
  // -- poissonic noise is only added if we have hits.
  if (hitList.size() > 0)
    for (XtalDiode xDiode; xDiode.isValid(); xDiode++) {
      DiodeNum diode = xDiode.getDiode();

      // convert dac in diode to MeV in xtal
      double meVXtal = (diode == LRG_DIODE) ?
        diodeDAC[xDiode]*mpdLrg.getVal() : // case LRG_DIODE:
        diodeDAC[xDiode]*mpdSm.getVal();   // case SM_DIODE:

      // MeV in xtal -> electrons in diode
      double eDiode = meVXtal * m_ePerMeV[diode];

      // apply poissonic fluctuation to # of electrons.
      double noise = sqrt(eDiode)*RandGauss::shoot();

      // add noise
      eDiode += noise;

      // convert back to dac in diode
      meVXtal = eDiode/m_ePerMeV[diode];
      diodeDAC[xDiode] = (diode == LRG_DIODE) ?
        meVXtal/mpdLrg.getVal() : // case LRG_DIODE:
        meVXtal/mpdSm.getVal();  // case SM_DIODE:
    }     // for xDiode

  
  // Store pedestal values here
  CalVec<XtalRng, float> pedVals(XtalRng::N_VALS);
  CalVec<XtalRng, float> pedSigs(XtalRng::N_VALS);

  // Stage 3: convert dac->adc
  for (XtalRng xRng; xRng.isValid(); xRng++) {
    XtalDiode xDiode = xRng.getXtalDiode();
    float cos;
    CalXtalId 
      rngXtalId(RngIdx(xtalId,xRng).getCalXtalId());
    

    // retrieve pedestals
    sc = m_calCalibSvc->getPed(rngXtalId,
                               pedVals[xRng],
                               pedSigs[xRng],
                               cos);
    if (sc.isFailure()) return sc;
    
    // calcuate adc val from dac
    sc = m_calCalibSvc->evalADC(rngXtalId, 
                                diodeDAC[xDiode], 
                                tmpADC[xRng]);
    if (sc.isFailure()) return sc;   
  } // for xRng

  // STAGE 4: electronic noise
  // electronic noise is calculated on a per-diode basis
  // uses sigmas from ADC pedestal data
  // adc values need not include pedestals for this calculation
  for (XtalDiode xDiode; xDiode.isValid(); xDiode++) {
    DiodeNum diode = xDiode.getDiode();

    // only add electronic nose to sm diodes if we have
    // hits, otherwise may cause false signals?
    if (diode == SM_DIODE && !hitList.size()) continue;

    FaceNum face = xDiode.getFace();
    
    float sigX8 = pedSigs[XtalRng(face,diode.getX8Rng())];
    float sigX1 = pedSigs[XtalRng(face,diode.getX1Rng())];
    
    // use same rand for both since the X1 & X8 noise
    // is strongly correlated
    double rnd = RandGauss::shoot();

    tmpADC[XtalRng(face,diode.getX1Rng())] += sigX1*rnd;
    tmpADC[XtalRng(face,diode.getX8Rng())] += sigX8*rnd;    
  }

  //-- LAC TESTS --//
  // operate on ped_subtraced adc
  CalibData::ValSig fle,fhe,lacThreshP,lacThreshN;
  // generate tmp xtalid's w/ POS & neg FACEs specified
  CalXtalId tmpIdP(xtalId.getTower(),
                   xtalId.getLayer(),
                   xtalId.getColumn(),
                   POS_FACE);
  CalXtalId tmpIdN(xtalId.getTower(),
                   xtalId.getLayer(),
                   xtalId.getColumn(),
                   NEG_FACE);
  
  // retreive thresholds for both faces
  sc = m_calCalibSvc->getTholdCI(tmpIdP,fle,fhe,lacThreshP);
  if (sc.isFailure()) return sc;
  //cout << "xtalId=" << (int)xtalId << "\ttmpIdP=" << (int)tmpIdP << "\ttmpIdN=" << (int)tmpIdN << endl;
  //cout << "xtalId=" << xtalId << "\ttmpIdP=" << tmpIdP << "\ttmpIdN=" << tmpIdN << endl;
  sc = m_calCalibSvc->getTholdCI(tmpIdN,fle,fhe,lacThreshN);
  if (sc.isFailure()) return sc;

  // set log-accept flags
  lacN = tmpADC[XtalRng(NEG_FACE,LEX8)] >= lacThreshN.getVal();
  lacP = tmpADC[XtalRng(NEG_FACE,LEX8)] >= lacThreshP.getVal();
  
  if (hitList.size() > 0) {
    lacN = tmpADC[XtalRng(NEG_FACE,LEX8)] >= lacThreshN.getVal();
    lacP = tmpADC[XtalRng(NEG_FACE,LEX8)] >= lacThreshP.getVal();
  }

  //-- ADD PEDS, HARDWARE RANGE CONSTRAINTS --//
  // double check that ADC vals are all >= 0
  // double check that ADC vals are all < maxadc(4095)
  for (XtalRng xRng; xRng.isValid(); xRng++) {
    tmpADC[xRng] += pedVals[xRng];
    tmpADC[xRng] = max<double>(0,tmpADC[xRng]);
    tmpADC[xRng] = min<double>(m_maxAdc,tmpADC[xRng]);
  }
    
  //-- BEST RANGE SELECTION --//
  // operates on adc + ped
  // pegged will be set true if HEX1 ADC saturates
  peggedP = peggedN = false; 
  for (FaceNum face; face.isValid(); face++) {
    // BEST RANGE
    RngNum rng;
    for (rng=0; rng.isValid(); rng++) {
      XtalRng xRng(face,rng);

      // get ULD threshold
      CalXtalId tmpId(xtalId.getTower(),
                      xtalId.getLayer(),
                      xtalId.getColumn(),
                      face, rng);
      CalibData::ValSig ULDThold;
      sc = m_calCalibSvc->getULDCI(tmpId,ULDThold);
      if (sc.isFailure()) return sc;

      // case of HEX1 range 
      // adc values have a ceiling in HEX1
      if (rng == RngNum::N_VALS-1) {
        if (tmpADC[xRng] > ULDThold.getVal()) {
          if (face == POS_FACE) peggedP = true;
          else peggedN = true;
          tmpADC[xRng] =  ULDThold.getVal();
          break; // quit before rng is incremented 1 too many
        }
      } else if (tmpADC[xRng] < ULDThold.getVal()) break;  // break on 1st energy rng that is < threshold
    }
    
    // assign range selection
    if (face == POS_FACE) rngP = (CalXtalId::AdcRange)(int)rng;
    else rngN = (CalXtalId::AdcRange)(int)rng;
  }
  
  //-- COPY ADCS TO OUTPUT --//
  // round to nearest integer.
  for (RngNum rng=0; rng.isValid(); rng++) {
    adcP[rng] = (int)floor(tmpADC[XtalRng(POS_FACE,rng)]+0.5);
    adcN[rng] = (int)floor(tmpADC[XtalRng(NEG_FACE,rng)]+0.5);
  }

  return StatusCode::SUCCESS;
}
