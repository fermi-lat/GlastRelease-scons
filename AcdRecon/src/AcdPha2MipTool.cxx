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
 { declareInterface<AcdIPha2MipTool>(this) ; }

AcdPha2MipTool::~AcdPha2MipTool()
{} 

// This function extracts geometry constants from xml
// file using GlastDetSvc
StatusCode AcdPha2MipTool::initialize()
 {
  MsgStream log(msgSvc(), name());
  StatusCode sc = StatusCode::SUCCESS;
  log<<MSG::INFO<<"BEGIN initialize()"<<endreq ;
 
  return sc;
 }

StatusCode AcdPha2MipTool::makeAcdHits (const Event::AcdDigiCol* digis,
					Event::AcdHitCol* hits)
  //
  //
  // TDS input:
  // TDS output:
{
  MsgStream log(msgSvc(),name()) ;
  output = hits;

  // loop on digis
  for ( Event::AcdDigiCol::const_iterator it = digis->begin();
	it != digis->end(); it++ ) {
    const Event::AcdDigi* aDigi = *it;
    if ( aDigi == 0 ) return StatusCode::FAILURE ;
    // do the intersections for this track

    float mipsPmtA(0.);
    float mipsPmtB(0.);
    
    bool calibOk = getCalibratedValues(*aDigi,mipsPmtA,mipsPmtB);
    if ( !calibOk ) return StatusCode::FAILURE;

    Event::AcdHit* newHit = new Event::AcdHit(*aDigi,mipsPmtA,mipsPmtB);
    if ( newHit == 0 ) return StatusCode::FAILURE;
    output->push_back(newHit);

  }
  // done
  return StatusCode::SUCCESS ;
}

bool AcdPha2MipTool::getCalibratedValues(const Event::AcdDigi& digi, float& mipsPmtA, float& mipsPmtB) const {

  // FIXME -- get these from calib DB
  static const float pedestal = 0.;  
  static const float mip = 1000.;       // this is in terms of counts above pedestal
  static const unsigned short FullScale = 0xFFFF;

  // do PMT A
  bool hasHitA = digi.getAcceptMapBit(Event::AcdDigi::A);
  Event::AcdDigi::Range rangeA = digi.getRange(Event::AcdDigi::A);  
  unsigned short phaA = hasHitA ? ( rangeA == Event::AcdDigi::LOW ? digi.getPulseHeight(Event::AcdDigi::A) : FullScale ) : 0.;
  float pedSubtracted_A = phaA -  pedestal;
  mipsPmtA = pedSubtracted_A / mip;

  // do PMT B
  bool hasHitB = digi.getAcceptMapBit(Event::AcdDigi::B);
  Event::AcdDigi::Range rangeB = digi.getRange(Event::AcdDigi::B);  
  unsigned short phaB = hasHitB ? ( rangeB == Event::AcdDigi::LOW ? digi.getPulseHeight(Event::AcdDigi::B) : FullScale ) : 0.;
  float pedSubtracted_B = phaB -  pedestal;
  mipsPmtB = pedSubtracted_B / mip;
 
  return true;
}
