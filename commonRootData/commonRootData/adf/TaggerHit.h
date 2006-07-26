
#ifndef ROOT_TaggerHit_H
#define ROOT_TaggerHit_H

#include "TObject.h"

/** @class TaggerHit
 * @brief The digitization ancillary data for beamtest 2006  
 * 
 * $Header$
*/

namespace commonRootData {

class TaggerHit: public TObject {

public:
    
    TaggerHit();

    TaggerHit(UInt_t moduleId, UInt_t layerId, UInt_t stripId,
             Double_t pulseHeight, Double_t sigma, Bool_t isPedSubtracted);

    TaggerHit(const TaggerHit& copy);

    void initialize(UInt_t moduleId, UInt_t layerId, UInt_t stripId,
             Double_t pulseHeight, Double_t sigma, Bool_t isPedSubtracted);
    
    TaggerHit& operator=(const TaggerHit& copy);

    virtual ~TaggerHit();

    void Clear(Option_t *option ="");

    void Print(Option_t *option="") const;

    UInt_t getModuleId() const { return m_moduleId; }
    UInt_t getLayerId() const { return m_layerId; }
    UInt_t getStripId() const { return m_stripId; }
    Double_t getPulseHeight() const { return m_pulseHeight; }
    Bool_t isPedestalSubtracted() const { return m_isPedestalSubtracted; } 
    Double_t getSigma() const { return m_sigma; }

    Bool_t CompareInRange( const TaggerHit &ref, const std::string& name="" )const;

private:

    UInt_t m_moduleId;
    UInt_t m_layerId;
    UInt_t m_stripId;
    Double32_t m_pulseHeight;
    Double32_t m_sigma;
    Bool_t m_isPedestalSubtracted;


    ClassDef(TaggerHit,1) // Digitization Ancillary data beamtest 2006
};

} // end namespace

#endif
