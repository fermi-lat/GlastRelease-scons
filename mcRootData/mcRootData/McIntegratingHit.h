#ifndef MCINTEGRATINGHIT_H
#define MCINTEGRATINGHIT_H

#include "TObject.h"
#include "TVector.h"
#include "TRef.h"

#include "VolumeIdentifier.h"
#include "McParticle.h"


#ifdef WIN32
using namespace std;
#endif

#include <map>

/** @class CompareTRef
 * @brief Defines comparision operator to use for TRef key in a map.
 *
 * @author Heather Kelly
 */
class CompareTRef {

public:
    CompareTRef() { };
    ~CompareTRef() { };

    bool operator() (const TRef& x, const TRef& y) const {
        return (x.Hash() > y.Hash());

    }
};

/** @class McIntegratingHit
* @brief
*
*  @author Heather Kelly
*  
* $Header$
*/
class McIntegratingHit: public TObject {
    
public:
    //typedef map<McParticle*,Double_t> energyDepositMap;
    typedef map<TRef, Double_t, CompareTRef> energyDepositMap;
    
    McIntegratingHit();
    
    virtual ~McIntegratingHit();

    void  Clear(Option_t *option ="");

    void initialize(VolumeIdentifier id, Double_t e, const TVector3& moment1, const TVector3& moment2);
    
    const VolumeIdentifier getVolumeId() const { return m_volumeId; };
    
    Double_t getTotalEnergy() const { return m_totalEnergy; };
    
    /// Retrieve itemized energy
    const energyDepositMap& getItemizedEnergy() const;
    /// Add single energyInfo to energyDepositMap
    void addEnergyItem( const double& energy, McParticle* t );
    
    /// Retrieve the energy-weighted first moments of the position
    const TVector3 getMoment1 () const;
    /// Retrieve the energy-weighted second moments of the position
    const TVector3 getMoment2 () const;
    
private:
        /// Packed flags for particle property
    UInt_t m_packedFlags;
    /// total deposited energy: set automatically when m_energyInfo is modified.
    Double_t m_totalEnergy;
    /// identifies what volume this integrating hit occurred in
    VolumeIdentifier m_volumeId;
    /// Vector of Energy information that consists of deposited energy and the mother McParticle
    energyDepositMap m_energyItem;
    /// Energy-weighted_first_moments_of_the_position * number_of_energy_deposition
    TVector3 m_moment1Seed;
    /// Energy-weighted_second_moments_of_the_position * number_of_energy_deposition
    TVector3 m_moment2Seed;
    
    ClassDef(McIntegratingHit,1)  // Monte Carlo Integrating Hit Class
};

#endif