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

  return sc;
 }

StatusCode AcdPha2MipTool::makeAcdHits (const Event::AcdDigiCol* digis,
					Event::AcdHitCol* hits)
  //
  //
  // TDS input:
  // TDS output:
{

  // sanity check
  if ( digis == 0 || hits == 0 ) {
    return StatusCode::FAILURE ;
  }

  // loop on digis  
  for ( Event::AcdDigiCol::const_iterator it = digis->begin();
	it != digis->end(); it++ ) {
    const Event::AcdDigi* aDigi = *it;
    if ( aDigi == 0 ) return StatusCode::FAILURE ;

    // get the calibrated values
    float mipsPmtA(0.);
    float mipsPmtB(0.);
    
    bool calibOk = getCalibratedValues(*aDigi,mipsPmtA,mipsPmtB);
    if ( !calibOk ) return StatusCode::FAILURE;

    Event::AcdHit* newHit = new Event::AcdHit(*aDigi,mipsPmtA,mipsPmtB);
    if ( newHit == 0 ) return StatusCode::FAILURE;

    hits->push_back(newHit);

  }
  // done
  return StatusCode::SUCCESS ;
}

bool AcdPha2MipTool::getCalibratedValues(const Event::AcdDigi& digi, float& mipsPmtA, float& mipsPmtB) const {

  static const unsigned short FullScale = 4095;

  // get calibration consts
  float pedValA(0.), pedValB(0.);
  float mipValA(0.), mipValB(0.);
				 
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
