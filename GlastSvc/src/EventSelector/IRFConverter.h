#ifndef IRFConverter_H 
#define IRFConverter_H 1

#include "instrument/DetectorConverter.h"
#include "instrument/Scintillator.h"
#include "instrument/SiDetector.h"
#include "instrument/CsIDetector.h"
#include "instrument/MCTruth.h"


#include "GlastEvent/TopLevel/ObjectVector.h"
#include "GlastEvent/MonteCarlo/MCACDhit.h"
#include "GlastEvent/MonteCarlo/MCTrack.h"
#include "GlastEvent/MonteCarlo/MCCalorimeterHit.h"

/* ! 
Derived from DetectorConverter, this class
provides a forward method to handle reading in of GlastDetector data 
Currently implmented is the ACD (Scintillator) into an ACDhitVector.  
*/

class IRFConverter : public DetectorConverter {
public:
	IRFConverter() {
		MCACDcontainer = new MCACDHitVector;
        MCCalorimeterHitContainer = new MCCalorimeterHitVector;
	}

    ~IRFConverter() { 
        if (MCACDcontainer)  delete MCACDcontainer;
        if (MCCalorimeterHitContainer) delete MCCalorimeterHitContainer;
    }

    // called due to a GlastDetector::accept(DetectorConverter) call
    virtual void forward (const Scintillator& s) {
        // ACD tile data
        // retrieve data only if this detector has data
        if ( !(s.empty()) ) {
            MCACDHit* tile = new MCACDHit(); 
            tile->setEnergy(s.energy());
            tile->setId(s.id());
            MCACDcontainer->push_back(tile);
        } 
    }

    virtual void forward ( const CsIDetector& csi) {
        // CAL CsI log data
        if ( !(csi.empty()) ) {
            
			MCCalorimeterHit* mcCal = new MCCalorimeterHit();
			mcCal->setEnergy(csi.energy());
			mcCal->setLayer(csi.layer());
			MCCalorimeterHitContainer->push_back(mcCal);

			// many more parameters need to be added here
        }
    }
    
    virtual void forward ( const MCTruth& mc) {
        // monte carlo truth
		/** Going to increase the functionality of
		    so that this converter can handle MCTrack class*/

        if ( !(mc.empty()) ) {

			/* Put on back burner untill Toby looks at it tomorrow

			MCTrack* particle = new MCTrack();
			HepLorentzVector* lv = new HepLorentzVector();
			lv->setX(mc.particle()->getIanCrap());
			lv->setY(mc.particle().m_pos.y());
			lv->setZ(mc.particle().m_pos.z());
			// need to get the set up the different 
			//properties here. Just for 

			particle->setFourMomentum(lv);
			MCTrackVector->push_back(particle); */

        }
    }

    virtual void forward ( const SiDetector& si) {
        // TKR silicon strip detectors
        if ( !(si.empty()) ) {

        }

    }

 // provide access to the list of ACD tile data
    MCACDHitVector* getACDTiles() { return MCACDcontainer; }
	MCCalorimeterHitVector* getMCCalHits() { return MCCalorimeterHitContainer; }

private:
    // one of Gaudi's ObjectContainers
	MCACDHitVector* MCACDcontainer;
	MCTrackVector* MCTrackContainer;
	MCCalorimeterHitVector* MCCalorimeterHitContainer;
    
};

#endif