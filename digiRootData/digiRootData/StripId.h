#ifndef StripId_H
#define StripId_H

#include "TObject.h"

/*! \class StripId
\brief Stores the data for a single hit TKR strip.  Provides all 
 necessary data to identify the location of this strip and the
 controller used to read it out - determines which ToT value this
 strip corresponds to.

 \li Jul 2001 Heather Kelly - provide storage for cable number as a 
          double-check for controller (left or right)
 \li Jun 2001 Heather Kelly - renamed StripID to StripId
 \li Jan 1999 Daniel Flath - ROOT HTML Documentation added
 \li Dec 1999 Daniel Flath - Creation
*/

class StripId : public TObject  
{
private:
    /*
    StripId tag word layout
    (defined in following enum):
    \verbatim
     _______________________________________________
    |15|  |  |12|11|10|  |  |  |  |  |  |  |  |  |00|
    |__|__|__|__|__|__|__|__|__|__|__|__|__|__|__|__|
    |   Tower   |CT|        Strip Number            |
    |___________|__|________________________________|
    \endverbatim

    ***** NOTES: *****

    (1) CT == Controller (right(1) or left (0))
    (2) CTRL, TOWER fields not used in '99/2000 beamtest - all were right controllers
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
    //! packed word containing strip data
    UShort_t m_tag;
    //! stores cable number [0-7]
    UChar_t m_cable; 
public:
    enum {
        STRIP_CTRL_LEFT = 0,
        STRIP_CTRL_RIGHT = 1
    };

    StripId();
    virtual ~StripId() { };

    //! returns the whole packed tag for this strip
    UShort_t getTag()   const { return m_tag; };
    //! returns strip id in [0,1599]
    UShort_t getId()    const { return (m_tag >> TKR_V_STRIP) & TKR_M_STRIP; };
    //! returns cable number used to readout this strip
    UChar_t  getCable() const { return m_cable; };
    //! returns the controller bit
    UShort_t getController()  const { return (m_tag >> TKR_V_CTRL) & TKR_M_CTRL; };
    //! returns the tower number
    UShort_t getTower() const { return (m_tag >> TKR_V_TOWER) & TKR_M_TOWER; };

    Bool_t setTag(UShort_t tagVal) { m_tag = tagVal; return kTRUE; };
    Bool_t setStrip(UShort_t stripVal);
    Bool_t setCable(UChar_t cableVal, UChar_t cntrl=0);
    Bool_t setController(UShort_t ctrlVal);
    Bool_t setTower(UShort_t towerVal);

    Int_t Compare(TObject *obj);
    Bool_t IsSortable() const;

    ClassDef(StripId,4)     // Information on a single Si Strip
};

#endif // !defined StripId_H
