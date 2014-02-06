#include "../AcdUtil/AcdUtilFuncs.h"


namespace AcdUtil {

  /**
   * @brief Fill an AcdEventToplogy object from a list of AcdHits


   * @param acdHits the hit
   * @param topology the topology object
   **/
  bool UtilFunctions::fillAcdEventTopology(const std::vector<Event::AcdHit*>& acdHits,
					   Event::AcdEventTopology& evtTopo){
    
    unsigned tileCount(0),ribbonCount(0),vetoCount(0),tileVeto(0);
    float totalTileEnergy(0),totalRibbonEnergy(0);  
    float tileEnergy(0),ribbonEnergy(0);  
    float ghostTileEnergy(0),ghostRibbonEnergy(0);  
    float triggerTileEnergy(0),triggerRibbonEnergy(0);  
    unsigned nTilesTop(0);  
    unsigned nTilesSideRow[4] = {0,0,0,0};  
    unsigned nTilesSideFace[4] = {0,0,0,0};  
    unsigned nVetoTop(0);  
    unsigned nVetoSideRow[4] = {0,0,0,0};  
    unsigned nVetoSideFace[4] = {0,0,0,0};  
    float tileEnergyTop(0);  
    float tileEnergySideRow[4] = {0.,0.,0.,0.};    
    float tileEnergySideFace[4] = {0.,0.,0.,0.};  
    unsigned nSidesHit(0),nSidesVeto(0);
    
    std::set<int> sidesHit;
    std::set<int> sidesVeto;
    
    for ( Event::AcdHitCol::const_iterator itr = acdHits.begin(); itr != acdHits.end(); itr++ ) {
      Event::AcdHit* theHit = *itr;
      const idents::AcdId& id = theHit->getAcdId();
      bool hasGhost   = theHit->getGhost();
      bool hasTrigger = theHit->getTriggerVeto();
      if ( hasTrigger)   vetoCount++;
      if ( id.tile() ) {
	tileCount++;
	sidesHit.insert( id.face() );
	if ( hasTrigger ) {
	  sidesVeto.insert( id.face() );
	  tileVeto++;
	  switch ( id.face() ) {
	  case 0:
	    nVetoTop++;
	    break;
	  case 1: case 2: case 3: case 4:
	    nVetoSideFace[id.face()-1]++;
	    nVetoSideRow[id.row()]++;
	  }
	}
	float energy = theHit->tileEnergy();
	totalTileEnergy += energy;
	tileEnergy += energy*!hasGhost;
	ghostTileEnergy += energy*hasGhost;
	triggerTileEnergy += energy*hasTrigger;
	switch ( id.face() ) {
	case 0:
	  nTilesTop++;
	  tileEnergyTop += energy;
	  break;
	case 1:
	case 2:
	case 3:
	case 4:
	  nTilesSideFace[id.face() - 1]++;
	  tileEnergySideFace[id.face() - 1] += energy;
	  nTilesSideRow[id.row()]++;
	  tileEnergySideRow[id.row()] += energy;
	  break;
	}      
      } else if ( id.ribbon() ) {
	ribbonCount++;
	float energy = theHit->ribbonEnergy( Event::AcdHit::A ) + theHit->ribbonEnergy( Event::AcdHit::B );
	energy /= 2.;
	totalRibbonEnergy += energy;
	ribbonEnergy += energy*!hasGhost;
	ghostRibbonEnergy += energy*hasGhost;
	triggerRibbonEnergy += energy*hasTrigger;
      }
    }
    
    nSidesHit = sidesHit.size();
    nSidesVeto = sidesVeto.size();
    evtTopo.set( tileCount,  ribbonCount,  vetoCount, tileVeto,
		 totalTileEnergy, totalRibbonEnergy, tileEnergy, ribbonEnergy,
               ghostTileEnergy, ghostRibbonEnergy, triggerTileEnergy, triggerRibbonEnergy,
		 nTilesTop,  nTilesSideRow,  nTilesSideFace,
		 nVetoTop,  nVetoSideRow,  nVetoSideFace,
		 tileEnergyTop,  tileEnergySideRow,  tileEnergySideFace,
		 nSidesHit,  nSidesVeto);
    
    return true;

  }

  /**
   * @brief fill cone energies from a set of AcdTkrHitPoca
   *
   * @param hitPocae the AcdTkrHitPoca
   * @param energy15 (energy in tiles within 15 degrees of the track)
   * @param energy30 (energy in tiles within 30 degrees of the track)
   * @param energy45 (energy in tiles within 45 degrees of the track)
   * @param triggerEnergy15 (energy within 15 degrees of the track for tiles with veto asserted)
   * @param triggerEnergy30 (energy within 30 degrees of the track for tiles with veto asserted)
   * @param triggerEnergy45 (energy within 45 degrees of the track for tiles with veto asserted)

   **/
  bool UtilFunctions::getConeEnergies(const std::vector<Event::AcdTkrHitPoca*>& hitPocae, 
				      float& energy15, float& energy30, float& energy45,
				      float& triggerEnergy15, float& triggerEnergy30, float& triggerEnergy45) {
    
    energy15 = energy30 = energy45 = 0.;
    triggerEnergy15 = triggerEnergy30 = triggerEnergy45 = 0.;
    
    static const double tan45 = 1.;
    static const double tan30 = 5.77350288616910401e-01;
    static const double tan15 = 2.67949200239410490e-01;
    
    for ( std::vector<Event::AcdTkrHitPoca*>::const_iterator itr = hitPocae.begin(); itr != hitPocae.end(); itr++ ) {
      const Event::AcdTkrHitPoca* hitPoca = *itr;
      // Don't include ribbons
      if ( hitPoca->getId().ribbon() ) continue;
      // Don't include stuff without hits
      if ( ! hitPoca->hasHit() ) continue;
      
      float tanAngle = ( -1 * hitPoca->getDoca() )/ hitPoca->getArcLength();
      float energy = hitPoca->tileEnergy();
      if (tanAngle < tan45) {
	if ( ! hitPoca->getGhost() ) {
	  energy45 += energy;
	}
	if ( hitPoca->getTriggerVeto() ) {
	  triggerEnergy45 += energy;
	}
	if (tanAngle < tan30) {
	  if ( ! hitPoca->getGhost() ) {
	    energy30 += energy;
	  }
	  if ( hitPoca->getTriggerVeto() ) {
	    triggerEnergy30 += energy;
	  }
	  if (tanAngle < tan15) {
	    if ( ! hitPoca->getGhost() ) {
	      energy15 += energy;
	    }
	    if ( hitPoca->getTriggerVeto() ) {
	      triggerEnergy15 += energy;
	  }
	  }
	}
      }
    }
  }


}
