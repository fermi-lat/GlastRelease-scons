
#ifndef overlayRootData_AcdOverlay_H
#define overlayRootData_AcdOverlay_H

#include "commonRootData/idents/AcdId.h"
#include "TObject.h"
#include "TString.h"
#include "TVector3.h"
#include "commonRootData/idents/VolumeIdentifier.h"

/** @class AcdOverlay
 * @brief The digitization data for a single ACD entity (tile or ribbon).  
 * 
 * $Header$
*/

class AcdOverlay: public TObject 
{
public:
    AcdOverlay();

    AcdOverlay(const VolumeIdentifier& volId);
    
    AcdOverlay(const VolumeIdentifier& volId, const AcdId& id, Double_t energy, const TVector3& position);
    
    virtual ~AcdOverlay() { };

    void Clear(Option_t *option ="");

    void Print(Option_t *option="") const;

    /// VolId for the detector
    const VolumeIdentifier& getVolId() const { return m_volumeId; };

    /// Returns the proper AcdId corresponding to the official numbering scheme
    const AcdId& getAcdId() const { return m_acdId; };

    /// Returns the estimate of energy deposited in this tile/ribbon
    Double_t getEnergyDep() const { return m_energyDep; };

    /// Returns the position
    const TVector3& getPosition() const {return m_center;}

    /// Returns the status word
    UInt_t getStatus() const {return m_status;}

    /// Root >= 3.0 is now const correct for the Compare function
    Int_t Compare(const TObject *obj) const; 
    Bool_t IsSortable() const;

private:

    /// Volume id for geometry
    VolumeIdentifier m_volumeId;
    /// ACD Id - contains methods to retrieve proper ACD id
    AcdId            m_acdId; 
    /// Energy deposited this hit
    Double_t         m_energyDep;
    /// Coordinate associated with this hit
    TVector3         m_center;
    /// status word
    UInt_t           m_status;

    ClassDef(AcdOverlay,1) // Digitization for a single ACD entity
};

#endif
