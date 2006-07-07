
#ifndef ROOT_QdcHit_H
#define ROOT_QdcHit_H

#include "TObject.h"

/** @class QdcHit
 * @brief The digitization ancillary data for beamtest 2006  
 * 
 * $Header$
*/

class QdcHit: public TObject {

public:
    
    QdcHit();

    QdcHit(UInt_t channel, UInt_t pulseHeight, Bool_t isPedSubtracted);

    QdcHit(const QdcHit& copy);

    void initialize(UInt_t channel, UInt_t pulseHeight, Bool_t isPedSubtracted);

    QdcHit& operator=(const QdcHit& copy);

    virtual ~QdcHit();

    void Clear(Option_t *option ="");

    void Print(Option_t *option="") const;

    UInt_t getChannel() const { return m_channel; }
    UInt_t getPulseHeight() const { return m_pulseHeight; }
    Bool_t isPedestalSubtracted() const { return m_isPedestalSubtracted; }

    Bool_t CompareInRange( const QdcHit &ref, const std::string& name="" ) const;

private:

    UInt_t m_channel;
    UInt_t m_pulseHeight;
    Bool_t m_isPedestalSubtracted;

    ClassDef(QdcHit,1) // Digitization Ancillary data beamtest 2006
};

#endif
