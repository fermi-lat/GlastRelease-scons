#ifndef IRFConverter_H 
#define IRFConverter_H 1

#include "instrument/DetectorConverter.h"

#include "idents/ModuleId.h"

//forward declarations

class TdCsIData;
#include "idents/ModuleId.h"
#include "GlastEvent/TopLevel/ObjectVector.h"
#include "GlastEvent/MonteCarlo/MCACDHit.h"
#include "GlastEvent/MonteCarlo/MCTrack.h"
#include "GlastEvent/MonteCarlo/MCCalorimeterHit.h"
#include "GlastEvent/MonteCarlo/MCSiLayer.h"
/*! 
Derived from DetectorConverter, this class
provides a forward method to handle reading in of GlastDetector data 
Currently implemented for :
    - MCACDHit (Scintillator) into an ACDhitVector.
    - MCCalorimeterHit (CsI Logs) into MCCalorimeterHitVector
    - MCTRKHit (Si Strips) int MCTRKHitVector
*/

class IRFConverter : public DetectorConverter {
public:
    
    //! constructor - create all object containers
    IRFConverter() ;

    //! destructor - delete all object containers
    ~IRFConverter();

    //! called due to a GlastDetector::accept(DetectorConverter) call
    virtual void forward (const Scintillator& s) ;

    //! called due to GlastDetector::accept(), handles Cal data
    virtual void forward ( const CsIDetector& csi) ;
    
    //! not implemented
    virtual void forward ( const MCTruth& mc);

    //! called due to GlastDetector::accept(), handles TKR strip data
    virtual void forward ( const SiDetector& si) ;

    //! sets m_tower
    virtual void forward( const Tower& t);

    //! provide access to the list of ACD tile data
    MCACDHitVector* getACDTiles() ;
    //! provide access to the list of CAL crystal data
    MCCalorimeterHitVector* getMCCalHits();
    //! provide access to the TKR strip data
    MCSiLayerVector* getMCTKRHits() ;
    //! provide access to the Raw CsI  Data geometry included as well
    TdCsIData* getTdCsIData() ;


private:
    idents::ModuleId m_towerId;

    // one of Gaudi's ObjectContainers
        TdCsIData* allcsiData;

	MCACDHitVector* MCACDcontainer;
	MCTrackVector* MCTrackContainer;
	MCCalorimeterHitVector* MCCalorimeterHitContainer;
        MCSiLayerVector* MCSiLayerContainer;  
    
};

#endif
