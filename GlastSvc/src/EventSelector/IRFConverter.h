#ifndef IRFConverter_H 
#define IRFConverter_H 1

#include "instrument/DetectorConverter.h"

#include "idents/ModuleId.h"

//forward declarations


#include "idents/ModuleId.h"
#include "GaudiKernel/ObjectVector.h"

#include "src/data/TdSiData.h"
#include "src/data/TdCsIData.h"
#include "src/data/TdVetoData.h"


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
    IRFConverter();
    
    //! destructor - delete all object containers
    ~IRFConverter(){};
    
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
    
    // these access functions for consistency with old
    const CsIData*  getCsIData()const{return &m_cal;}
    const SiData*   getSiData()const{return &m_tkr;}
    const IVetoData* getVetoData()const{return &m_acd;}
    
    void clear();
    
    //! print out contents
    virtual void printOn(std::ostream& out)const;
   
    
private:
    idents::ModuleId m_id;
    
    TdCsIData  m_cal;
    TdSiData   m_tkr;
    TdVetoData m_acd;
    
};

#endif
