#ifndef overlayRootData_TkrOverlay_H
#define overlayRootData_TkrOverlay_H 1

#include "TObject.h"
#include <vector>
#ifndef R__GLOBALSTL
#ifndef WIN32
using std::vector;
#else
using namespace std;
#endif
#endif

#include "commonRootData/idents/TowerId.h"
#include "GlastAxis.h"

/**
* @class TkrOverlay
* @brief PDS version of the TKR overlay objects
*
* Contains layer and tower identification, ToT and a list of hit strips
* The digis are produced either from MC hit output or from the actual data
* --> essentially a clone of the TkrDigi objects
* @author Leon Rochester
* $Header$
*/

class TkrOverlay : public TObject 
{    
public:
    
    TkrOverlay();
    
    ~TkrOverlay ();
    
    void initialize (Int_t l, GlastAxis::axis v, TowerId t, Int_t* tot);
    
    void Clear(Option_t *option ="");

    void Print(Option_t *option="") const;
    
    /**@defgroup TkrOverlayGroup TkrOverlay End-User Interface */
    /*@{*/
    //! Retrieve tower info
    TowerId getTower() const { return m_tower; }
    //! Retrieve layer info
    UShort_t getBilayer() const { return m_bilayer; } 
    //! Retrieve view info
    GlastAxis::axis getView() const { return m_view; }
    //! Retrieve the last strip associated with controller zero
    Int_t getLastController0Strip() const { return m_lastController0Strip; }
    //! Retrieve ToT info
    Int_t getToT(UInt_t i) const;
    //! Get ToT for a given strip
    Int_t getToTForStrip(UInt_t strip) const;
    //! number of hits
    UInt_t getNumHits() const { return m_hitCol.size(); }
    //! Get the ith strip number
    Int_t getStrip(UInt_t i) const;
    //! Get the ith strip
    Int_t getHit(UInt_t i) const;

    //! Get the pointer to the hit list
    const vector<UInt_t> getHitCol() const { return m_hitCol; };
    /*@}*/

    //! Add a controller 0 hit to the hit list (keeps track of highest strip)
    void addC0Hit( UInt_t strip );
    
    //! Add a controller 1 hit to the hit list
    void addC1Hit( UInt_t strip);

//    void* operator new(size_t size);

//    void* operator new(size_t size, void* vp);

//    void  operator delete(void* p);
    
private:

    UShort_t m_bilayer;
    UShort_t m_tot[2];
    /// Denotes highest strip number from controller 0, set to -1 when undefined.
    Short_t m_lastController0Strip; 
    TowerId m_tower;
    /// Denotes which axis is measured
    GlastAxis::axis m_view;
    /// Collection of hit strips
    vector<UInt_t> m_hitCol;
    
    ClassDef(TkrOverlay,1)
};

#endif // overlayRootData_TkrOverlay_H
