#ifndef ACDUTILFUNCS_H
#define ACDUTILFUNCS_H

#include "Event/Recon/AcdRecon/AcdHit.h"
#include "Event/Recon/AcdRecon/AcdTkrHitPoca.h"
#include "Event/Recon/AcdRecon/AcdEventTopology.h"

namespace AcdUtil {

  namespace UtilFunctions {

    bool fillAcdEventTopology(const std::vector<Event::AcdHit*>& acdHits,
			      Event::AcdEventTopology& topology);

    
    bool getConeEnergies(const std::vector<Event::AcdTkrHitPoca*>& hitPocae, 
			 float& energy15, float& energy30, float& energy45,
			 float& triggerEnergy15, float& triggerEnergy30, float& triggerEnergy45);
    

  }
}

#endif
