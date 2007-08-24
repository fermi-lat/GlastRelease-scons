#include "AcdPha2MipTool.h"
#include "GaudiKernel/MsgStream.h"

#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/DeclareFactoryEntries.h"

#include "Event/Recon/AcdRecon/AcdHit.h"
#include "Event/Digi/AcdDigi.h"

#include "idents/AcdId.h"

DECLARE_TOOL_FACTORY(AcdPha2MipTool)

AcdPha2MipTool::AcdPha2MipTool
 ( const std::string & type, 
   const std::string & name,
   const IInterface * parent )
 : AlgTool( type, name, parent )
 { 
   declareInterface<AcdIPha2MipTool>(this) ; 
   declareProperty("AcdCalibSvc",    m_calibSvcName = "AcdCalibSvc");
   declareProperty("PHATileCut",    m_pha_tile_cut = 0.0);
   declareProperty("MIPSTileCut",    m_mips_tile_cut = 0.0);
   declareProperty("PHARibbonCut",    m_pha_ribbon_cut = 0.0);
   declareProperty("MIPSRibbonCut",    m_mips_ribbon_cut = 0.0);
 }

AcdPha2MipTool::~AcdPha2MipTool()
{} 

// This function extracts geometry constants from xml
// file using GlastDetSvc
StatusCode AcdPha2MipTool::initialize()
 {
  MsgStream log(msgSvc(), name());
  StatusCode sc = StatusCode::SUCCESS;

  sc = service(m_calibSvcName, m_calibSvc, true);
  if ( !sc.isSuccess() ) {
    log << MSG::ERROR << "Could not get CalibDataSvc " << m_calibSvcName << endreq;
    return sc;
  } else {
    log << MSG::INFO << "Got CalibDataSvc " << m_calibSvcName << endreq;
  }

  log<<MSG::INFO<<"BEGIN initialize()"<<endreq ;
  
  return sc;
 }

StatusCode AcdPha2MipTool::makeAcdHits( const Event::AcdDigiCol& digis,
					bool periodicEvent, 
					Event::AcdHitCol& hits,
					AcdRecon::AcdHitMap& hitMap)
  //
  //
  // TDS input:
  // TDS output:
{


  MsgStream log(msgSvc(),name()) ;

  // loop on digis
  for ( Event::AcdDigiCol::const_iterator it = digis.begin();
	it != digis.end(); it++ ) {
    const Event::AcdDigi* aDigi = *it;
    // sanity check
    if ( aDigi == 0 ) return StatusCode::FAILURE ;
    
    // get the hit mask bits
    unsigned int hitMask = 0;
    hitMask |= aDigi->getAcceptMapBit(Event::AcdDigi::A) ? AcceptMapBit_AMask : 0;
    hitMask |= aDigi->getAcceptMapBit(Event::AcdDigi::B) ? AcceptMapBit_BMask : 0;
    hitMask |= aDigi->getVeto(Event::AcdDigi::A) ? VetoBit_AMask : 0;
    hitMask |= aDigi->getVeto(Event::AcdDigi::B) ? VetoBit_BMask : 0;
    hitMask |= aDigi->getCno(Event::AcdDigi::A) ? CNO_AMask : 0;
    hitMask |= aDigi->getCno(Event::AcdDigi::B) ? CNO_BMask : 0;	
    hitMap[aDigi->getId()] = hitMask;
    
    Event::AcdHit* newHit(0);
    StatusCode sc = makeAcdHit(*aDigi,periodicEvent,newHit);

    if ( sc.isFailure() ) return sc;
    if ( newHit != 0 ) {
      hits.push_back(newHit);
    }
  } 
  return StatusCode::SUCCESS ;
}

StatusCode AcdPha2MipTool::makeAcdHit ( const Event::AcdDigi& digi,
					bool periodicEvent, 
					Event::AcdHit*& hit) {
  float mipsPmtA(0.);
  float mipsPmtB(0.);
  bool acceptDigi(false);
  bool ok = getCalibratedValues(digi,mipsPmtA,mipsPmtB,acceptDigi);
  if ( !ok ) return StatusCode::FAILURE;
  if ( acceptDigi ) {
    hit = new Event::AcdHit(digi,mipsPmtA,mipsPmtB);
  } else {
    hit = 0;
  }
  return StatusCode::SUCCESS;
}

bool acceptAcdDigi ( const Event::AcdDigi& digi, bool periodicEvent );


bool AcdPha2MipTool::getCalibratedValues(const Event::AcdDigi& digi, float& mipsPmtA, float& mipsPmtB, bool& acceptDigi) const {

  static const unsigned short FullScale = 4095;

  // get calibration consts
  float pedValA(0.), pedValB(0.);
  float mipValA(0.), mipValB(0.);
				 
  acceptDigi = false;

  if ( ! getPeds(digi.getId(),pedValA,pedValB) ) {
    MsgStream log(msgSvc(), name());
    log << MSG::ERROR << "Failed to get pedestals." << endreq;
    return false;
  }
  if ( ! getMips(digi.getId(),mipValA,mipValB) ) {
    MsgStream log(msgSvc(), name());
    log << MSG::ERROR << "Failed to get gains." << endreq;
    return false;
  }

  // do PMT A
  bool hasHitA = digi.getAcceptMapBit(Event::AcdDigi::A);
  Event::AcdDigi::Range rangeA = digi.getRange(Event::AcdDigi::A);  
  unsigned short phaA = hasHitA ? ( rangeA == Event::AcdDigi::LOW ? digi.getPulseHeight(Event::AcdDigi::A) : FullScale ) : 0;
  if ( phaA == 0 ) {
    mipsPmtA = 0.;
  } else {
    float pedSubtracted_A = phaA -  pedValA;    
    mipsPmtA = pedSubtracted_A / mipValA;
    acceptDigi |= accept(digi.getId(),pedSubtracted_A,mipsPmtA);
  }

  // do PMT B
  bool hasHitB = digi.getAcceptMapBit(Event::AcdDigi::B);
  Event::AcdDigi::Range rangeB = digi.getRange(Event::AcdDigi::B);  
  unsigned short phaB = hasHitB ? ( rangeB == Event::AcdDigi::LOW ? digi.getPulseHeight(Event::AcdDigi::B) : FullScale ) : 0;
  if ( phaB == 0 ) {
    mipsPmtB = 0.;
  } else {
    float pedSubtracted_B = phaB -  pedValB;
    mipsPmtB = pedSubtracted_B / mipValB;
    acceptDigi |= accept(digi.getId(),pedSubtracted_B,mipsPmtB);
  } 

  return true;
}


bool AcdPha2MipTool::getPeds(const idents::AcdId& id, float& valA, float& valB) const {
  if ( m_calibSvc == 0 ) return false;  
  CalibData::AcdPed* ped(0);

  StatusCode sc = m_calibSvc->getPedestal(id,Event::AcdDigi::A,ped);
  if ( sc.isFailure() ) {
    return false;
  }
  valA = ped->getMean();

  sc = m_calibSvc->getPedestal(id,Event::AcdDigi::B,ped);
  if ( sc.isFailure() ) {
    return false;
  }
  valB = ped->getMean();
  
  return true;
}



bool AcdPha2MipTool::getMips(const idents::AcdId& id, float& valA, float& valB) const {
  if ( m_calibSvc == 0 ) return false;

  CalibData::AcdGain* gain(0);

  StatusCode sc = m_calibSvc->getMipPeak(id,Event::AcdDigi::A,gain);
  if ( sc.isFailure() ) {
    return false;
  }
  valA = gain->getPeak();
  
  sc = m_calibSvc->getMipPeak(id,Event::AcdDigi::B,gain);
  if ( sc.isFailure() ) {
    return false;
  }
  valB = gain->getPeak();
  
  return true;
}

bool AcdPha2MipTool::accept(const idents::AcdId& id, float pedSubtracted, float mips) const {
  if ( id.tile() ) {
    if ( pedSubtracted < m_pha_tile_cut ) return false;
    if ( mips < m_mips_tile_cut ) return false;    
  } else if ( id.ribbon() ) {
    if ( pedSubtracted < m_pha_ribbon_cut ) return false;
    if ( mips < m_mips_ribbon_cut ) return false;        
  } else {
    return false;
  }
  return true;
}
