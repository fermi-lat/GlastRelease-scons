
#ifndef ROOT_ScalerHit_H
#define ROOT_ScalerHit_H

#include "TObject.h"
#include <string>

/** @class ScalerHit
 * @brief The digitization ancillary data for beamtest 2006  
 * 
 * $Header$
*/

namespace commonRootData {

class ScalerHit: public TObject {

public:
    
    ScalerHit();

    ScalerHit(UInt_t channel, UInt_t val);

    ScalerHit(const ScalerHit& copy);

    void initialize(UInt_t ch, UInt_t val);
    
    ScalerHit& operator=(const ScalerHit& copy);

    virtual ~ScalerHit();

    void Clear(Option_t *option ="");

    void Print(Option_t *option="") const;

    UInt_t getChannel() const { return m_channel; }
    UInt_t getValue() const { return m_value; }

    Bool_t CompareInRange( const ScalerHit &ref, const std::string& name="" )const;

private:

    UInt_t m_channel;
    UInt_t m_value;


    ClassDef(ScalerHit,1) // Digitization Ancillary data beamtest 2006
};

} // end namespace

#endif
