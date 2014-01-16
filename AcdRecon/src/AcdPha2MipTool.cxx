#include "AcdPha2MipTool.h"
#include "GaudiKernel/MsgStream.h"

#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/DeclareFactoryEntries.h"

#include "Event/Recon/AcdRecon/AcdHit.h"
#include "Event/Digi/AcdDigi.h"

#include "CalibData/Acd/AcdPed.h"
#include "CalibData/Acd/AcdGain.h"
#include "CalibData/Acd/AcdHighRange.h"
#include "CalibData/Acd/AcdCoherentNoise.h"

#include "AcdUtil/AcdCalibFuncs.h"

#include "idents/AcdId.h"

DECLARE_TOOL_FACTORY(AcdPha2MipTool)

AcdPha2MipTool::AcdPha2MipTool
 ( const std::string & type, 
   const std::string & name,
   const IInterface * parent )
   : AlgTool( type, name, parent ),
     m_gemDeltaEventTime(0)
 { 
   declareInterface<AcdIPha2MipTool>(this) ; 
   declareProperty("AcdCalibSvc",    m_calibSvcName = "AcdCalibSvc");
   declareProperty("PHATileCut",    m_pha_tile_cut = 45.0);  // Set to 45 to match flight data
   declareProperty("MIPSTileCut",    m_mips_tile_cut = 0.0);
   declareProperty("PHARibbonCut",    m_pha_ribbon_cut = 45.0);  // Set to 45 to match flight data
   declareProperty("MIPSRibbonCut",    m_mips_ribbon_cut = 0.0);
   declareProperty("VetoThrehsold",    m_vetoThreshold = 0.4);
   declareProperty("ApplyCoherentNoiseCalib",    m_applyCoherentNoiseCalib = false);
     
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
					unsigned gemDeltaEventTime, 
					Event::AcdHitCol& hits,
					AcdRecon::AcdHitMap& hitMap)
  //
  //
  // TDS input:
  // TDS output:
{


  MsgStream log(msgSvc(),name()) ;

  m_gemDeltaEventTime = gemDeltaEventTime;

  // loop on digis
  for ( Event::AcdDigiCol::const_iterator it = digis.begin();
	it != digis.end(); it++ ) {
    const Event::AcdDigi* aDigi = *it;
    // sanity check
    if ( aDigi == 0 ) return StatusCode::FAILURE ;
    
    // Don't do with NA channels, they are likely to have bad calibrations, which can cause problems
    if ( ! aDigi->getId().tile() && ! aDigi->getId().ribbon() ) continue;

    // get the hit mask bits
    //unsigned int hitMask = 0;
    //hitMask |= aDigi->getAcceptMapBit(Event::AcdDigi::A) ? AcceptMapBit_AMask : 0;
    //hitMask |= aDigi->getAcceptMapBit(Event::AcdDigi::B) ? AcceptMapBit_BMask : 0;
    //hitMask |= aDigi->getVeto(Event::AcdDigi::A) ? VetoBit_AMask : 0;
    //hitMask |= aDigi->getVeto(Event::AcdDigi::B) ? VetoBit_BMask : 0;
    //hitMask |= aDigi->getCno(Event::AcdDigi::A) ? CNO_AMask : 0;
    //hitMask |= aDigi->getCno(Event::AcdDigi::B) ? CNO_BMask : 0;	
    
    Event::AcdHit* newHit(0);
    StatusCode sc = makeAcdHit(*aDigi,periodicEvent,newHit);

    if ( sc.isFailure() ) return sc;
    if ( newHit != 0 ) {
      hits.push_back(newHit);
      hitMap[aDigi->getId()] = newHit;
    }

 
  } 
  return StatusCode::SUCCESS ;

}

StatusCode AcdPha2MipTool::makeAcdHit ( const Event::AcdDigi& digi,
					bool /* periodicEvent */, 
					Event::AcdHit*& hit) {
  double mipsPmtA(0.);
  double mipsPmtB(0.);
  bool acceptDigi(false);
  bool ok = getCalibratedValues(digi,mipsPmtA,mipsPmtB,acceptDigi);

  // Check Veto thresholds
  if ( digi.getHitMapBit( Event::AcdDigi::A ) ) {
    mipsPmtA = mipsPmtA > m_vetoThreshold ? mipsPmtA : m_vetoThreshold;
    acceptDigi = true;
  }
  if ( digi.getHitMapBit( Event::AcdDigi::B ) ) {
    mipsPmtB = mipsPmtB > m_vetoThreshold ? mipsPmtB : m_vetoThreshold;
    acceptDigi = true;
  }

  // Check for Ninja Hits
  if ( digi.isNinja() || digi.getGemFlag() ) {
    mipsPmtA = mipsPmtA > m_vetoThreshold ? mipsPmtA : m_vetoThreshold;
    mipsPmtB = mipsPmtB > m_vetoThreshold ? mipsPmtB : m_vetoThreshold;
    acceptDigi = true;
  }

  if ( !ok ) return StatusCode::FAILURE;
  if ( acceptDigi ) {
    hit = new Event::AcdHit(digi,mipsPmtA,mipsPmtB);
  } else {
    hit = 0;
  }
  return StatusCode::SUCCESS;
}


bool AcdPha2MipTool::getCalibratedValues(const Event::AcdDigi& digi, double& mipsPmtA, double& mipsPmtB, bool& acceptDigi) const {

  // get calibration consts
  acceptDigi = false;
  double pedSubA(0.);
  double pedSubB(0.);

  // do PMT A
  bool hasHitA = digi.getAcceptMapBit(Event::AcdDigi::A) || digi.getVeto(Event::AcdDigi::A);
  if ( hasHitA ) {
    Event::AcdDigi::Range rangeA = digi.getRange(Event::AcdDigi::A);  
    bool ok = rangeA == Event::AcdDigi::LOW ? 
      getValues_lowRange(digi.getId(),Event::AcdDigi::A,digi.getPulseHeight(Event::AcdDigi::A),pedSubA,mipsPmtA) :
      getValues_highRange(digi.getId(),Event::AcdDigi::A,digi.getPulseHeight(Event::AcdDigi::A),pedSubA,mipsPmtA);
    if ( !ok ) return false;
    acceptDigi |= rangeA == Event::AcdDigi::HIGH ? true : accept(digi.getId(),pedSubA,mipsPmtA);
  }

  // do PMT B
  bool hasHitB = digi.getAcceptMapBit(Event::AcdDigi::B) || digi.getVeto(Event::AcdDigi::A);
  if ( hasHitB ) {
    Event::AcdDigi::Range rangeB = digi.getRange(Event::AcdDigi::B);  
    bool ok = rangeB == Event::AcdDigi::LOW ? 
      getValues_lowRange(digi.getId(),Event::AcdDigi::B,digi.getPulseHeight(Event::AcdDigi::B),pedSubB,mipsPmtB) :
      getValues_highRange(digi.getId(),Event::AcdDigi::B,digi.getPulseHeight(Event::AcdDigi::B),pedSubB,mipsPmtB);
    if ( !ok ) return false;
    acceptDigi |= rangeB == Event::AcdDigi::HIGH ? true : accept(digi.getId(),pedSubB,mipsPmtB);
  }

  return true;
}

bool AcdPha2MipTool::getValues_lowRange(const idents::AcdId& id, Event::AcdDigi::PmtId pmt, unsigned short pha, 
					double& pedSub, double& mips) const {
  
  if ( m_calibSvc == 0 ) return false;  
  CalibData::AcdPed* ped(0);
  
  StatusCode sc = m_calibSvc->getPedestal(id,pmt,ped);
  if ( sc.isFailure() ) {
    return false;
  }
  double pedestal = ped->getMean();

  CalibData::AcdGain* gain(0);
  sc = m_calibSvc->getMipPeak(id,pmt,gain);
  if ( sc.isFailure() ) {
    return false;
  }
  double mipPeak = gain->getPeak();

  if ( m_applyCoherentNoiseCalib ) {
    CalibData::AcdCoherentNoise* cNoise(0);
    sc = m_calibSvc->getCoherentNoise(id,pmt,cNoise);
    if ( sc.isFailure() ) {
      return false;
    }
    double deltaPed(0.);
    sc = AcdCalib::coherentNoise(m_gemDeltaEventTime,
				 cNoise->getAmplitude(),cNoise->getDecay(),cNoise->getFrequency(),cNoise->getPhase(),
				 deltaPed);
    if ( sc.isFailure() ) {
      return false;
    }
    pedestal += deltaPed;
  }

  pedSub = (double)pha - pedestal;
  sc = AcdCalib::mipEquivalent_lowRange(pha,pedestal,mipPeak,mips);
  return sc.isFailure() ? false : true;

}

bool AcdPha2MipTool::getValues_highRange(const idents::AcdId& id, Event::AcdDigi::PmtId pmt, unsigned short pha, 
					 double& pedSub, double& mips) const {

  if ( m_calibSvc == 0 ) return false;  
  CalibData::AcdHighRange* highRange(0);
  StatusCode sc = m_calibSvc->getHighRange(id,pmt,highRange);
  if ( sc.isFailure() ) {
    return false;
  }
  double pedestal = highRange->getPedestal();
  double slope = highRange->getSlope();
  double saturation = highRange->getSaturation();
  pedSub = (double)pha - pedestal;
  sc = AcdCalib::mipEquivalent_highRange(pha,pedestal,slope,saturation,mips);
  return sc.isFailure() ? false : true;

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
