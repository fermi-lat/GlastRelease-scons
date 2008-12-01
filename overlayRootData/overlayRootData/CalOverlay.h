
#ifndef overlayRootData_CalOverlay_H
#define overlayRootData_CalOverlay_H 1

#include "commonRootData/idents/CalXtalId.h"
#include "TObject.h"
#include "TVector3.h"

/** @class CalOverlay
 * @brief Root CAL data object for overlays
 *
 * @author Tracy Usher
 * $Header$
*/

class CalOverlay : public TObject  { 
    
public:
    
    /** @defgroup CalOverlayGroup CalOverlay End-User Interface */
    /*@{*/
    /// Status bit definitions
    enum DigiStatusBits {
         DIGI_NULL    = 0,
         DIGI_MC      = 1<<28,   // This digi is from Monte Carlo
         DIGI_DATA    = 1<<30,   // This digi is from data 
         DIGI_OVERLAY = 1<<31    // This digi contains overlay info
    };
    /*@}*/
    
    CalOverlay();
    
    virtual ~CalOverlay();

    void initialize(const CalXtalId &id, const TVector3& position, const Double_t& energy, const UInt_t& status);

    void Clear(Option_t *option ="");

    void Print(Option_t *option="") const;

    /** @ingroup CalOverlayGroup */
    /*@{*/
    /// Retrieve log identifier
    const CalXtalId getPackedId() const { return m_xtalId;   };
    const TVector3& getPosition() const { return m_position; };
    const Double_t  getEnergy()   const { return m_energy;   };
    const UInt_t    getStatus()   const { return m_status;   };
    /*@}*/

    void addToStatus(DigiStatusBits& bits2Add) {m_status |= bits2Add;}

private:
    
    CalXtalId m_xtalId;
    TVector3  m_position;
    Double_t  m_energy;
    UInt_t    m_status;
    
    ClassDef(CalOverlay,1) // CAL Overlay Class
};

#endif
