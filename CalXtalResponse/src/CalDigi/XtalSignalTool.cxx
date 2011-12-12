//    $Header$

/** @file implement XtalSignalTool.h
    @author Z.Fewtrell
*/

// Include files

// LOCAL
#include "XtalSignalTool.h"
#include "CalXtalResponse/ICalCalibSvc.h"

// GLAST
#include "Event/MonteCarlo/McIntegratingHit.h"
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "idents/CalXtalId.h"
#include "CalibData/Cal/CalMevPerDac.h"
#include "idents/VolumeIdentifier.h"
#include "Event/TopLevel/Event.h"
#include "Event/TopLevel/EventModel.h"
#include "CalUtil/CalDefs.h"


// EXTLIB
#include "CLHEP/Random/RandGauss.h"

#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/IDataProviderSvc.h"

#include "TTree.h"

// STD
#include <map>
#include <string>

const float XtalSignalTool::m_ePerMeVInDiode = 2.77e5;

//static ToolFactory<XtalSignalTool> s_factory;
//const IToolFactory& XtalSignalToolFactory = s_factory;
DECLARE_TOOL_FACTORY(XtalSignalTool);

/// declare jobOptions properties
XtalSignalTool::XtalSignalTool( const std::string& type,
                                const std::string& name,
                                const IInterface* parent)
  : AlgTool(type,name,parent),
    m_calCalibSvc(0),
    m_nCsISeg(0),
    m_CsILength(0),
    m_eXtal(0),
    m_eDiodeMSm(0),
    m_eDiodePSm(0),
    m_eDiodeMLarge(0),
    m_eDiodePLarge(0)
{
  declareInterface<IXtalSignalTool>(this);

  declareProperty("CalCalibSvc",        m_calCalibSvcName    = "CalCalibSvc");
}

StatusCode XtalSignalTool::initialize() {
  MsgStream msglog(msgSvc(), name());   
  msglog << MSG::INFO << "initialize" << endreq;

  StatusCode sc;


  //-- jobOptions --//
  if ((sc = setProperties()).isFailure()) {
    msglog << MSG::ERROR << "Failed to set properties" << endreq;
    return sc;
  }

  //-- Retrieve constants from GlastDetSVc
  sc = retrieveConstants();
  if (sc.isFailure())
    return sc;

  ///////////////////////////////////////
  //-- RETRIEVE HELPER TOOLS & SVCS  --//
  ///////////////////////////////////////

  // obtain CalCalibSvc
  sc = service(m_calCalibSvcName.value(), m_calCalibSvc);
  if (sc.isFailure()) {
    msglog << MSG::ERROR << "can't get " << m_calCalibSvcName  << endreq;
    return sc;
  }

  return StatusCode::SUCCESS;
}


StatusCode XtalSignalTool::retrieveConstants() {
  double val;
  typedef std::map<int*,string> PARAMAP;

  // try to find the GlastDetSvc service
  IGlastDetSvc* detSvc;
  StatusCode sc = service("GlastDetSvc", detSvc);
  if (sc.isFailure() ) {
    MsgStream msglog(msgSvc(), name());   
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

  for(PARAMAP::iterator it=param.begin(); it!=param.end();it++)
    if(!detSvc->getNumericConstByName((*it).second, &val)) {
      MsgStream msglog(msgSvc(), name());   
      msglog << MSG::ERROR << " constant " <<(*it).second <<" not defiPed" << endreq;
      return StatusCode::FAILURE;
    } else *((*it).first)=(int)val;

  // doubles are done separately
  double tmp;
  sc = detSvc->getNumericConstByName("CsILength", &tmp);
  m_CsILength = tmp;
  if (sc.isFailure()) {
    MsgStream msglog(msgSvc(), name());   
    msglog << MSG::ERROR << " constant CsILength not defiPed" << endreq;
    return sc;
  }
  
  return StatusCode::SUCCESS;
}

/**
   convert energy deposit to diode signal levels in CIDAC (charge injection DAC) units.
   Use asymmetry and MeVPerDAC calibrations to calculate signal level for each diode.
 */
StatusCode XtalSignalTool::calculate(const Event::McIntegratingHit &hit,
                                     CalUtil::CalVec<CalUtil::XtalDiode, float> &cidacArray) {
  StatusCode sc;

  m_dat.Clear();

  const idents::CalXtalId xtalId((idents::VolumeIdentifier)hit.volumeID());
  m_dat.xtalIdx = CalUtil::XtalIdx(xtalId);

  /////////////////////////////////////
  // STAGE 1: fill signal energies ////
  /////////////////////////////////////

  //-- retrieve mevPerDAC calib
  CalibData::CalMevPerDac const * const mpd = m_calCalibSvc->getMPD(m_dat.xtalIdx);
  if (!mpd) return StatusCode::FAILURE;
  m_dat.mpd[CalUtil::LRG_DIODE] = mpd->getBig()->getVal();
  m_dat.mpd[CalUtil::SM_DIODE]  = mpd->getSmall()->getVal();
      
  // get volume identifier.
  const idents::VolumeIdentifier volId((idents::VolumeIdentifier)hit.volumeID());
      
  // make sure the hits are cal hits
  if (!volId.isCal()) {
     MsgStream msglog( msgSvc(), name() );
     msglog << MSG::WARNING << "mcHit is not Cal." << endreq ;
     return StatusCode::SUCCESS;
  }

    
  //--XTAL DEPOSIT--//
  if((int)volId[CalUtil::fCellCmp] ==  m_eXtal ) {
    sc = sumCsIHit(hit, cidacArray);
    if (sc.isFailure()) return sc;
  } 
  //--DIODE DEPOSIT--//
  else {
    sc = sumDiodeHit(hit, cidacArray);
    if (sc.isFailure()) return sc;
  }
  
  return StatusCode::SUCCESS;
}

/**
   Algorithm has following steps:
   - determine mean cidac signal using mevPerDAC calibrations
   - sum indivual diode signal levels in propotion to light asymmetry calibration
 */
StatusCode XtalSignalTool::sumCsIHit(const Event::McIntegratingHit &hit,
                                     CalUtil::CalVec<CalUtil::XtalDiode, float> &cidacArray) {

  StatusCode sc;

  const float ene = hit.totalEnergy();
  
  const float realpos = hit2pos(hit);

  //--ASYMMETRY--//
  float asymL, asymS;
      
  sc = m_calCalibSvc->evalAsym(m_dat.xtalIdx, CalUtil::ASYM_LL, realpos, asymL);
  if (sc.isFailure()) return sc;

  sc = m_calCalibSvc->evalAsym(m_dat.xtalIdx, CalUtil::ASYM_SS, realpos, asymS);
  if (sc.isFailure()) return sc;

  // Get CIDAC values for each diode
  const Event::McIntegratingHit::XtalEnergyDep& xtalDepDSN = hit.getXtalEnergyDep("diodeSNEG");
  double diodeSNTotE = xtalDepDSN.getTotalEnergy();
  double diodeSNDirE = xtalDepDSN.getDirectEnergy();
  //double cidacSN     = (diodeSNTotE + diodeSNDirE) / m_dat.mpd[CalUtil::SM_DIODE];
  double cidacSN     = diodeSNTotE / m_dat.mpd[CalUtil::SM_DIODE];
  double cidacSNDir  = diodeSNDirE / m_dat.mpd[CalUtil::SM_DIODE];

  const Event::McIntegratingHit::XtalEnergyDep& xtalDepDLN = hit.getXtalEnergyDep("diodeLNEG");
  double diodeLNTotE = xtalDepDLN.getTotalEnergy();
  double diodeLNDirE = xtalDepDLN.getDirectEnergy();
//  double cidacLN     = (diodeLNTotE + diodeLNDirE) / m_dat.mpd[CalUtil::LRG_DIODE];
  double cidacLN     = diodeLNTotE / m_dat.mpd[CalUtil::LRG_DIODE];
  double cidacLNDir    = diodeLNDirE / m_dat.mpd[CalUtil::LRG_DIODE];

  const Event::McIntegratingHit::XtalEnergyDep& xtalDepDSP = hit.getXtalEnergyDep("diodeSPOS");
  double diodeSPTotE = xtalDepDSP.getTotalEnergy();
  double diodeSPDirE = xtalDepDSP.getDirectEnergy();
//  double cidacSP     = (diodeSPTotE + diodeSPDirE) / m_dat.mpd[CalUtil::SM_DIODE];
  double cidacSP     = diodeSPTotE / m_dat.mpd[CalUtil::SM_DIODE];
  double cidacSPDir     = diodeSPDirE / m_dat.mpd[CalUtil::SM_DIODE];

  const Event::McIntegratingHit::XtalEnergyDep& xtalDepDLP = hit.getXtalEnergyDep("diodeLPOS");
  double diodeLPTotE = xtalDepDLP.getTotalEnergy();
  double diodeLPDirE = xtalDepDLP.getDirectEnergy();
//  double cidacLP     = (diodeLPTotE + diodeLPDirE) / m_dat.mpd[CalUtil::LRG_DIODE];
  double cidacLP     = diodeLPTotE / m_dat.mpd[CalUtil::LRG_DIODE];
  double cidacLPDir     = diodeLPDirE / m_dat.mpd[CalUtil::LRG_DIODE];

  //-- DETERMINE MEAN CIDAC SIGNAL LEVEL --//
//  const float meanCIDACS = ene/m_dat.mpd[CalUtil::SM_DIODE];
//  const float meanCIDACL = ene/m_dat.mpd[CalUtil::LRG_DIODE];

  //  calc cidacSmP, cidacSmN, cidacLrgP, cidacLrgN - here are quick notes
  //   asym=log(p/n)
  //   mean=sqrt(p*N)
  //   exp(asym)*mean^2=p/n*p*n=p^2
  //   p=exp(asym/2)*mean
  //   n=exp(-asym/2)*mean

  //-- SUM ENERGY TO EACH DIODE --//
  using namespace CalUtil;
//  cidacArray[XtalDiode(POS_FACE, LRG_DIODE)] += exp(   asymL/2) * meanCIDACL;
//  cidacArray[XtalDiode(NEG_FACE, LRG_DIODE)] += exp(-1*asymL/2) * meanCIDACL;
//  cidacArray[XtalDiode(POS_FACE, SM_DIODE)]  += exp(   asymS/2) * meanCIDACS;
//  cidacArray[XtalDiode(NEG_FACE, SM_DIODE)]  += exp(-1*asymS/2) * meanCIDACS;
  cidacArray[XtalDiode(POS_FACE, LRG_DIODE)] += exp(   asymL/2) * cidacLP + 1.5*cidacLPDir;
  cidacArray[XtalDiode(NEG_FACE, LRG_DIODE)] += exp(-1*asymL/2) * cidacLN + 1.5*cidacLNDir;
  cidacArray[XtalDiode(POS_FACE, SM_DIODE)]  += exp(   asymS/2) * cidacSP + 1.5*cidacSPDir;
  cidacArray[XtalDiode(NEG_FACE, SM_DIODE)]  += exp(-1*asymS/2) * cidacSN + 1.5*cidacSNDir;

#if 0 // debug printing
  cout << "XtalSignalTool::sumCsIHit(): " << m_dat.xtalIdx.toStr() 
	   << " " << ene
	   << " " << realpos
	   << " " << asymL
	   << " " << asymS
       << endl;
  for (XtalDiode xDiode;xDiode.isValid(); xDiode++)
	  cout << " " << xDiode.val() << " " << cidacArray[xDiode] << endl;
#endif

  return StatusCode::SUCCESS;
}

float XtalSignalTool::hit2pos(const Event::McIntegratingHit &hit) {
  // get volume identifier.
  const idents::VolumeIdentifier volId((idents::VolumeIdentifier)hit.volumeID());

  //-- HIT POSITION--//
  const HepPoint3D mom1 = hit.moment1();
  const int segm = volId[CalUtil::fSegment]; // segment # (0-11)D
  // let's define the position of the segment along the crystal
  float relpos = (segm+0.5)/m_nCsISeg; // units in xtal len 0.0-1.0
  // in local reference system x is always oriented along the crystal
  const float dpos =  mom1.x();

  // add position w/in segment to segment position
  relpos += dpos/m_CsILength; // still ranges 0.0-1.0

  // in 'mm' from center of xtal (units used in asym splines)
  return (relpos-0.5)*m_CsILength;
}

StatusCode XtalSignalTool::sumDiodeHit(const Event::McIntegratingHit &hit,
                                       CalUtil::CalVec<CalUtil::XtalDiode, float> &cidacArray) {

  StatusCode sc;
        
  float ene = hit.totalEnergy();
  float eDiode = ene*m_ePerMeVInDiode; // convert MeV-in-diode to electrons

  // get volume identifier.
  idents::VolumeIdentifier volId = ((idents::VolumeIdentifier)hit.volumeID());

  CalUtil::DiodeNum diode;
  CalUtil::FaceNum face;
  if((int)volId[CalUtil::fCellCmp]      == m_eDiodePLarge) face = CalUtil::POS_FACE, diode = CalUtil::LRG_DIODE;
  else if((int)volId[CalUtil::fCellCmp] == m_eDiodeMLarge) face = CalUtil::NEG_FACE, diode = CalUtil::LRG_DIODE;
  else if((int)volId[CalUtil::fCellCmp] == m_eDiodePSm)    face = CalUtil::POS_FACE, diode = CalUtil::SM_DIODE;
  else if((int)volId[CalUtil::fCellCmp] == m_eDiodeMSm)    face = CalUtil::NEG_FACE, diode = CalUtil::SM_DIODE;
  else {
      MsgStream msglog(msgSvc(), name());   
      msglog << MSG::ERROR << "VolId is not in xtal or diode.  Programmer error." << endreq;
      return StatusCode::FAILURE;
  }
      
  // convert e-in-diode to MeV-in-xtal-center
  const float meVXtal = eDiode/m_ePerMeV[diode.val()]; 

  // convert MeV deposition to CIDAC val.
  // use asymmetry to determine how much of mean cidac to apply to diode in question
  // remember we are treating mevXtal as at the xtal center (pos=0.0)
  float meanCIDAC;
  float asym;
  meanCIDAC = meVXtal/m_dat.mpd[diode];
  if (diode == CalUtil::LRG_DIODE)
    sc = m_calCalibSvc->getAsymCtr(m_dat.xtalIdx, CalUtil::ASYM_LL, asym);
  else
    sc = m_calCalibSvc->getAsymCtr(m_dat.xtalIdx, CalUtil::ASYM_SS, asym);
  if (sc.isFailure()) return sc;

  // meanCIDAC = sqrt(cidacP*cidacN), exp(asym)=P/N, P=N*exp(asym), N=P/exp(asym)
  // mean^2 = P*N, mean^2=P*P/exp(asym), mean^2*exp(asym)=P^2
  float cidacP;
  cidacP = sqrt(meanCIDAC*meanCIDAC*exp(asym));
  const float cidac = (face==CalUtil::POS_FACE) ? cidacP : cidacP/exp(asym);

  // sum cidac val to running total
  const CalUtil::XtalDiode xDiode(face,diode);
  cidacArray[xDiode] += cidac;

  return StatusCode::SUCCESS;
}

StatusCode XtalSignalTool::finalize() {

  return StatusCode::SUCCESS;
}
