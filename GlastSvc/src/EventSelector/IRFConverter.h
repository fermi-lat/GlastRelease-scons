#ifndef IRFConverter_H 
#define IRFConverter_H 1

#include "instrument/DetectorConverter.h"

#include "idents/ModuleId.h"

//forward declarations


#include "idents/ModuleId.h"
#include "GlastEvent/TopLevel/ObjectVector.h"

#include "GlastEvent/Irf/IrfAcdHit.h"
#include "GlastEvent/Irf/IrfCalHit.h"
#include "GlastEvent/Irf/IrfTkrLayer.h"
#include "src/Raw/LdSiData.h"
#include "src/Raw/LdCsIData.h"

/*! 
Derived from DetectorConverter, this class
provides a forward method to handle reading in of GlastDetector data 
Currently implemented for :
    - IrfAcdHit (Scintillator) into an IrfAcdHitVector.
    - IrfCalHit (CsI Logs) into IrfCalHitVector
    - IrfTkrHit (Si Strips) int IrfTkrHitVector
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

    //! provide access to the Raw CsI  Data geometry included as well
    TdCsIData* getTdCsIData() ;

    TdSiData* getTdSiData();

    //! provide access to the list of Acd hits
    IrfAcdHitVector* getIrfAcdHits();
    //! provide access to the list of Cal hits
    IrfCalHitVector* getIrfCalHits();
    //! provide access to the list of Tkr hits
    IrfTkrLayerVector* getIrfTkrHits();

private:
    idents::ModuleId m_towerId;

    // one of Gaudi's ObjectContainers
        LdCsIData* allcsiData;
        LdSiData* allsiData;


        IrfAcdHitVector* IrfAcdHitContainer;
        IrfCalHitVector* IrfCalHitContainer;
        IrfTkrLayerVector* IrfTkrLayerContainer;
    
};

#endif
