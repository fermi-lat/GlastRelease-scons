#ifndef MCINTEGRATINGHIT_H
#define MCINTEGRATINGHIT_H

#include "TObject.h"
#include "TVector.h"
#include "TRefArray.h"
#include "TMap.h"

#include "commonRootData/idents/VolumeIdentifier.h"
#include "McParticle.h"

#include <vector>
#ifndef R__GLOBALSTL
#ifndef WIN32
using std::vector;
#else
using namespace std;
#endif
#endif


/** @class McIntegratingHit
* @brief Stores a single McIntegratingHit.  This class mirrors the TDS
* Event::McIntegratingHit class.
*
* Each McIntegratingHit is composed of:
* - VolumeIdentifier
* - Status flags packed into an unsigned integer
* - Total energy
* - First and Second moments
* - One of the following:
*    -# An array of energies for the primary particle (and e-/e+ if available)
*    -# A map of McParticle pointers and their corresponding energies
*
* If the map of McParticleIds and energy are available, they are stored as 2 
* separate arrays.  One is a vector of Int_t for the ids, and the other a 3
* element array of energy for the primary particle, and optionally if the 
* primary is gamma which pair produces, an e-/e+ pair.  Access is provided
* through methods that retrieve the appropriate McParticleId, energy pair.
* Use the methods itemizedEnergyIdSize(), itemizedEnergyIdReset(), and
* Int_t itemizedEnergyIdNext(Double_t& energy) to access the data.
*
* If available, the map of McParticles and energies are actually stored to the 
* ROOT file as 2 separate arrays.  One TRefArray of McParticles and a vector of
* energies.  Access is provided through methods that retrieve the appropriate
* McParticle*, energy pair.  Use the methods itemizedEnergySize(), 
* itemizedEnergyReset() and 
* const* McParticle itemizedEnergyNext(Double_t& energy) to access the data.
*
* @author Heather Kelly
*  
* $Header$
*/


/** @page mcRootData_doc
    @section McIntegratingHit McIntegratingHit Documentation
    Notes

<table>
<tr><th><EM> Method <th> Return Type </EM><th> Description <th>

<small>
<tr><td> getVolumeId
<td>const VolumeIdentifier<td> Returns the volume identifier
<tr><td> getTotalEnergy
<td>Double_t<td> Returns the total energy deposited in a volume
<tr><td> getMoment1
<td>const TVector3<td> Returns the first moment 
<tr><td> getMoment2
<td>const TVector3<td> Returns the second moment
<tr><td> getMcParticleEnergy( Particle p )
<td>Double_t<td> Returns energy associated with PRIMARY, ELECTRON or POSITRON
</small>
</table>
*/


// Define sub class to contain energy information for each diode when in a crystal volume
class McXtalEnergyDep : public TObject
{
public:
    McXtalEnergyDep() : m_totalEnergy(0.), m_directEnergy(0.), m_moment1seed(0.,0.,0.), m_moment2seed(0.,0.,0.) {}
   ~McXtalEnergyDep() {}

    // Apparently, we need a copy constructor to keep compiler happy
    McXtalEnergyDep(const McXtalEnergyDep& dep) {*this = dep;}

    void Clear(Option_t * ="");
    void Fake( Int_t ievent, UInt_t rank, Float_t randNum ) ; // for tests
    Bool_t CompareInRange( const McXtalEnergyDep&, const std::string& name = "" ) const ; // for tests
    void Print(Option_t * ="") const;

    /// Get total deposited energy seen by this diode
    Double_t getTotalEnergy()  const {return m_totalEnergy;}

    /// Get the direct deposited energy seen by this diode
    Double_t getDirectEnergy() const {return m_directEnergy;}

    /// Retrieve the energy-weighted first moments of the position
    const TVector3 getMoment1() const;
          TVector3 getMoment1();

    /// Retrieve the energy-weighted second moments of the position
    const TVector3 getMoment2() const;
          TVector3 getMoment2();

    // add energy/position for this diode
    void addEnergyItem(Double_t totalEnergy, Double_t directEnergy, const TVector3& position);

    // For initialization
    void initialize(Double_t totalEnergy, Double_t directEnergy, const TVector3& mom1, const TVector3& mom2);

private:
    /// total deposited energy: set automatically when m_energyInfo is modified.
    Double_t  m_totalEnergy;
    /// "direct" energy seen by diode
    Double_t  m_directEnergy;
    /// Energy-weighted_first_moments_of_the_position * number_of_energy_deposition
    TVector3  m_moment1seed;
    /// Energy-weighted_second_moments_of_the_position * number_of_energy_deposition
    TVector3  m_moment2seed;

    ClassDef(McXtalEnergyDep,1)  // XtalEnergyDep Class
};

class McIntegratingHit: public TObject {
    
public:
    
    /// This enum is used to aid in the fill of the m_energyArray member
    enum Particle{
        PRIMARY = 0,
        ELECTRON = 1,
        POSITRON = 2
    };

    McIntegratingHit();
    
    virtual ~McIntegratingHit();

    // Apparently, we need a copy constructor to keep compiler happy
    McIntegratingHit(const McIntegratingHit& hit) {*this = hit;}

    void Clear(Option_t * ="");
    void Fake( Int_t ievent, UInt_t rank, Float_t randNum ) ; // for tests
    Bool_t CompareInRange( const McIntegratingHit &, const std::string & name = "" ) const ; // for tests
    void Print(Option_t * ="") const;

    void initialize(const VolumeIdentifier& id);
    
    const VolumeIdentifier getVolumeId() const { return m_volumeId; };
    
    Double_t getTotalEnergy() const { return m_totalEnergy; };
    
    /// Add single energyInfo to McParticle*, energy map
    void addEnergyItem( const Double_t& energy, McParticle* t, const TVector3& pos );

    void setEnergyItems( const Double_t &totE, const Double_t *energyArr, 
        const TVector3& moment1, const TVector3& moment2);

    /// Retrieve the next McParticle, energy pair
    /// Returns the pair corresponding to m_mapPtr
    /// If we have finished traversing the list - Next will return null (0).
    const McParticle* itemizedEnergyNext (Double_t &energy) const ;
    /// Set the McParticle, energy map counter to zero - the beginning
    void itemizedEnergyReset() const { m_mapPtr = 0; };
    /// Returns the size of the list of McParticle, energy pairs
    UInt_t itemizedEnergySize() const { return m_mcPartArr.GetEntries(); };

    /// RetriEve energy associated with PRIMARY, ELECTRON or POSITRON
    Double_t getMcParticleEnergy(Particle p) const;

    /// Retrieve the energy-weighted first moments of the position
    const TVector3 getMoment1 () const;
    /// Retrieve the energy-weighted second moments of the position
    const TVector3 getMoment2 () const;

    /// Retrieve XtalEnergyDep object
    const McXtalEnergyDep* getXtalEnergyDep(const TString& key) const;
    McXtalEnergyDep* getXtalEnergyDep(const TString& key);

    /// Direct access to map for filling
    TMap& getXtalEnergyDepMap() {return m_xtalEnergyMap;}

    void* operator new(size_t size);

    void* operator new(size_t size, void* vp);

    McIntegratingHit& operator=(const McIntegratingHit& rhs);
    
private:
    /// total deposited energy: set automatically when m_energyInfo is modified.
    Double_t         m_totalEnergy;
    /// Optionally store the energy associated with the primary and e-/e+
    Double_t         m_energyArray[3];
    /// Packed flags for particle property
    UInt_t           m_packedFlags;
    /// identifies what volume this integrating hit occurred in
    VolumeIdentifier m_volumeId;
    /// Energy-weighted_first_moments_of_the_position * number_of_energy_deposition
    TVector3         m_moment1Seed;
    /// Energy-weighted_second_moments_of_the_position * number_of_energy_deposition
    TVector3         m_moment2Seed;
    
    /// Actually store McParticles as a TRefArray
    TRefArray        m_mcPartArr;
    /// Store energies from McParticle* pair in a vector
    vector<Double_t> m_energyPtrArr;

    /// Mapping between diodes and deposited energy
    TMap             m_xtalEnergyMap;
   
    /// Store a vector of McParticle Ids
    //std::vector<Int_t> m_mcIdArr;

    /// Vector of Energy information that consists of deposited energy and 
    /// the mother McParticle.  This data member is transient - not written!
    //energyDepositMap m_energyItem; //!
    
    /// transient data member to keep track of place in McParticle*, energy map
    mutable UInt_t m_mapPtr; //!

    ClassDef(McIntegratingHit,2)  // Monte Carlo Integrating Hit Class
};

#endif
