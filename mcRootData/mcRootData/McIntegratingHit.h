#ifndef MCINTEGRATINGHIT_H
#define MCINTEGRATINGHIT_H

#include "TObject.h"
#include "TVector.h"
#include "TRefArray.h"

#include "VolumeIdentifier.h"
#include "McParticle.h"

#ifdef WIN32
using namespace std;
#endif

#include <vector>

/** @class McIntegratingHit
* @brief Stores a single McIntegratingHit.  This class mirrors the TDS
* McIntegratingHit class.
* Each McIntegratingHit is composed of:
* - VolumeIdentifier
* - Status flags packed into an unsigned integer
* - Total energy
* - First and Second moments
* - A map of McParticles and their corresponding energies
*
* The map of McParticles and energies are actually stored to the ROOT file
* as 2 separate arrays.  One TRefArray of McParticles and a vector of
* energies.  Access is provided through methods that retrieve the appropriate
* McParticle*, energy pair.  Use the methods mapSize(), mapReset() and 
* const* McParticle mapNext(Double_t& energy) to acces the data.
*
* @author Heather Kelly
*  
* $Header$
*/
class McIntegratingHit: public TObject {
    
public:
    
    McIntegratingHit();
    
    virtual ~McIntegratingHit();

    void Clear(Option_t *option ="");

    void Print(Option_t *option="") const;

    void initialize(const VolumeIdentifier& id);
    
    const VolumeIdentifier getVolumeId() const { return m_volumeId; };
    
    Double_t getTotalEnergy() const { return m_totalEnergy; };
    
    /// Add single energyInfo to energyDepositMap
    void addEnergyItem( const Double_t& energy, McParticle* t, const TVector3& pos );
    /// Retrieve the next McParticle, energy pair
    /// Returns the pair corresponding to m_mapEntry - set via the mapReset()
    /// If we have finished traversing the list - mapNext will return null (0).
    const McParticle* mapNext(Double_t &energy);
    /// Set the McParticle, energy map counter to zero - the beginning
    void mapReset() { m_mapEntry = 0; };
    /// Returns the size of the list of McParticle, energy pairs
    UInt_t mapSize() const { return m_mcPartArr.GetEntries(); };

    /// Retrieve the energy-weighted first moments of the position
    const TVector3 getMoment1 () const;
    /// Retrieve the energy-weighted second moments of the position
    const TVector3 getMoment2 () const;
    
private:
    /// total deposited energy: set automatically when m_energyInfo is modified.
    Double_t m_totalEnergy;
    /// Packed flags for particle property
    UInt_t m_packedFlags;
    /// identifies what volume this integrating hit occurred in
    VolumeIdentifier m_volumeId;
    /// Energy-weighted_first_moments_of_the_position * number_of_energy_deposition
    TVector3 m_moment1Seed;
    /// Energy-weighted_second_moments_of_the_position * number_of_energy_deposition
    TVector3 m_moment2Seed;
    /// Actually store McParticles as a TRefArray
    TRefArray m_mcPartArr;
    /// Store energies in a vector
    std::vector<Double_t> m_energies;
    /// Vector of Energy information that consists of deposited energy and 
    /// the mother McParticle.  This data member is transient - not written!
    //energyDepositMap m_energyItem; //!
    UInt_t m_mapEntry; //!

    ClassDef(McIntegratingHit,1)  // Monte Carlo Integrating Hit Class
};

#endif