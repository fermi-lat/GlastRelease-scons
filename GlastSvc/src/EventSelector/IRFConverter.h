#ifndef IRFConverter_H 
#define IRFConverter_H 1

#include "instrument/DetectorConverter.h"
#include "instrument/Scintillator.h"
#include "instrument/SiDetector.h"
#include "instrument/CsIDetector.h"

#include "GlastEvent/TopLevel/ObjectVector.h"
#include "GlastEvent/Hits/ACDhit.h"

/* ! 
Derived from DetectorConverter, this class
provides a forward method to handle reading in of GlastDetector data 
Currently implmented is the ACD (Scintillator) into an ACDhitVector.  
*/

class IRFConverter : public DetectorConverter {
public:
	IRFConverter() {
		ACDcontainer = new ACDhitVector;
	}

    ~IRFConverter() { if (ACDhitVector)  delete ACDcontainer;}

    // called due to a GlastDetector::accept(DetectorConverter) call
    virtual void forward (const Scintillator& s) {
        // ACD tile data
        // retrieve data only if this detector has data
        if ( !(s.empty()) ) {
            ACDhit* tile = new ACDhit();
            tile->setEnergy(s.energy());
            tile->setId(s.id());
            ACDcontainer->push_back(tile);
        } 
    }

    virtual void forward ( const CsIDetector& csi) {
        // CAL CsI log data
        if ( !(csi.empty()) ) {

        }
    }
    
    virtual void forward ( const MCTruth& mc) {
        // monte carlo truth
        if ( !(mc.empty()) ) {

        }
    }

    virtual void forward ( const SiDetector& si) {
        // TKR silicon strip detectors
        if ( !(si.empty()) ) {

        }

    }

 // provide access to the list of ACD tile data
    ACDhitVector* getACDTiles() { return ACDcontainer; }

private:
    // one of Gaudi's ObjectContainers
	ACDhitVector* ACDcontainer;
    
};

#endif