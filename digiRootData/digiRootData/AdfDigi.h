
#ifndef ROOT_AdfDigi_H
#define ROOT_AdfDigi_H

#include "TObject.h"
#include "TClonesArray.h"

#include "TaggerHit.h"
#include "QdcHit.h"

/** @class AdfDigi
 * @brief The digitization ancillary data for beamtest 2006  
 * 
 * $Header$
*/

class AdfDigi: public TObject {

public:
    
    AdfDigi();

    AdfDigi(const AdfDigi &copy);

    AdfDigi& operator=(const AdfDigi& copy);

    virtual ~AdfDigi();

    void Clear(Option_t *option ="");

    void Print(Option_t *option="") const;

    void setEventNumber(UInt_t evtNum) { m_eventNumber = evtNum;}
    void setSpillNumber(UInt_t spillNum) { m_spillNumber = spillNum; }

    UInt_t getEventNumber() const { return m_eventNumber; }
    UInt_t getSpillNumber() const { return m_spillNumber; }

    const TClonesArray* getTaggerHitCol() const { return m_taggerHitCol; }
    const TaggerHit* getTaggerHit(UInt_t ind) const;
    TaggerHit* addTaggerHit(UInt_t moduleId, UInt_t layerId, UInt_t stripId,
                            UInt_t pulseHgt, Bool_t isPedSubtracted);


    const TClonesArray* getQdcHitCol() const { return m_qdcHitCol; }
    const QdcHit* getQdcHit(UInt_t ind) const;
    QdcHit* addQdcHit(UInt_t channel, UInt_t pulseHgt, Bool_t isPedSub);


    void Fake( Int_t ievent, Float_t randNum );
    Bool_t CompareInRange( AdfDigi &, const std::string & name = "" );

private:

    UInt_t m_eventNumber;
    UInt_t m_spillNumber;

    Int_t m_numTaggerHit;
    Int_t m_numQdcHit;

    TClonesArray *m_taggerHitCol;
    static TClonesArray *s_taggerHitStaticCol; //!

    TClonesArray *m_qdcHitCol;
    static TClonesArray *s_qdcHitStaticCol; //!

    ClassDef(AdfDigi,1) // Digitization Ancillary data beamtest 2006
};

#endif
