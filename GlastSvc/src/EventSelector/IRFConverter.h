#ifndef IRFConverter_H 
#define IRFConverter_H 1

#include "instrument/DetectorConverter.h"
#include "instrument/Scintillator.h"

#include "GlastEvent/TopLevel/ObjectVector.h"
#include "GlastEvent/Hits/ACDhit.h"

/* ! 
Derived from DetectorConverter, this class
provides a forward method to handle reading in of ACD data (Scintillator)
into an ACDhitVector.  
*/

class IRFConverter : public DetectorConverter {
public:
	IRFConverter() {
		ACDcontainer = new ACDhitVector;
	}

    ~IRFConverter() { if (ACDhitVector)  delete ACDcontainer;}

    // called due to a GlastDetector::accept(DetectorConverter) call
    virtual void forward (const Scintillator& s) {
        if ( !(s.empty()) ) {
            ACDhit* tile = new ACDhit();
            tile->setEnergy(s.energy());
            tile->setId(s.id());
            ACDcontainer->push_back(tile);
        } 
    }

 // provide access to the list of ACD tile data
    ACDhitVector* getACDTiles() { return ACDcontainer; }

private:
	ACDhitVector* ACDcontainer;
    
};

#endif