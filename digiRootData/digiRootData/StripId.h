#ifndef StripId_H
#define StripId_H

#include "TObject.h"

/*! \class StripId
\brief
 Jun 2001 Heather Kelly - renamed StripID to StripId
 Jan 1999 Daniel Flath - ROOT HTML Documentation added
 Dec 1999 Daniel Flath - Creation
*/

class StripId : public TObject  
{
private:
    /*
    StripId tag word layout
    (defined in following enum):
     _______________________________________________
    |15|  |  |12|11|10|  |  |  |  |  |  |  |  |  |00|
    |__|__|__|__|__|__|__|__|__|__|__|__|__|__|__|__|
    |   Tower   |CT|        Strip Number            |
    |___________|__|________________________________|

    ***** NOTES: *****

    (1) CT == Controller (right or left)
    (2) CTRL, TOWER fields not used in '99/2000 beamtest
    (3) Valid strip numbers: [0,1599] for full layers
    (4) CTRL field is used in the 2001 balloon, still only one tower

    */
    enum {
        TKR_K_STRIP = 11,
        TKR_V_STRIP = 0,
        TKR_M_STRIP = ((1 << TKR_K_STRIP) - 1),
        TKR_K_CTRL = 1,
        TKR_V_CTRL = (TKR_K_STRIP + TKR_V_STRIP),
        TKR_M_CTRL = ((1 << TKR_K_CTRL) - 1),
        TKR_K_TOWER = 4,
        TKR_V_TOWER = (TKR_K_CTRL + TKR_V_CTRL),
        TKR_M_TOWER = ((1 << TKR_K_TOWER) - 1)
    };
    UShort_t m_tag;     // Packed word containing strip data
public:
    enum {
        STRIP_CTRL_RIGHT = 0,
        STRIP_CTRL_LEFT = 1
    };

    StripId();
    virtual ~StripId() { };

    UShort_t getTag()   const { return m_tag; };
    UShort_t getId()    const { return (m_tag >> TKR_V_STRIP) & TKR_M_STRIP; };
    UShort_t getController()  const { return (m_tag >> TKR_V_CTRL) & TKR_M_CTRL; };
    UShort_t getTower() const { return (m_tag >> TKR_V_TOWER) & TKR_M_TOWER; };

    Bool_t setTag(UShort_t tagVal) { m_tag = tagVal; return kTRUE; };
    Bool_t setStrip(UShort_t stripVal);
    Bool_t setController(UShort_t ctrlVal);
    Bool_t setTower(UShort_t towerVal);

    Int_t Compare(TObject *obj);
    Bool_t IsSortable() const;

    ClassDef(StripId,3)     // Information on a single Si Strip
};

#endif // !defined StripId_H
