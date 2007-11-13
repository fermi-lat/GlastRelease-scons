#ifndef MCTRAJECTORY_H
#define MCTRAJECTORY_H

#include "TObject.h"
#include "TObjArray.h"
#include "TClonesArray.h"

#include "commonRootData/idents/VolumeIdentifier.h"
#include "McParticle.h"

/** @class McTrajectoryPoint
* @brief GLAST Monte Carlo Particle Trajectory "point" class.
*
* The McTrajectoryPoint class contains the step energy, volume Id and position information 
* for each step taken by a particle being tracked through GLAST during the simulation stage.
* - A volume identifier to give the volume this step is in
* - The energy of the track at this step
* - A TVector with the position
*  
*/


/** @page mcRootData_doc
    @section McTrajectoryPoint McTrajectoryPoint Documentation
    Notes

<table>
<tr><th><EM> Method <th> Return Type </EM><th> Description <th>

<small>
<tr><td> getVolumeID
<td>VolumeIdentifier<td> Returns the volume identifier
<tr><td> getEnergy
<td>Float_t<td> Returns the energy at this  point 
<tr><td> getPoint
<td>TVector3<td> Returns the hit point 
</small>
</table>
*/

class McTrajectoryPoint : public TObject
{
public:
    McTrajectoryPoint();
    McTrajectoryPoint(VolumeIdentifier vId, float energy, TVector3& point) :
        m_volumeID(vId), m_energy(energy), m_point(point) {}
    ~McTrajectoryPoint() {}

    /// clear lists, free pointers, etc., after read from / write to file
    void Clear(Option_t *option ="");
    void Fake( Int_t ievent, UInt_t rank, Float_t randNum ) ; // for tests
    Bool_t CompareInRange( const McTrajectoryPoint &, const std::string & name = "" ) const ; // for tests
    void Print(Option_t *option="") const;

    /// initialize 
    void initialize(VolumeIdentifier id, Float_t energy, const TVector3& position);

    /// Retrieve cell identifier
    VolumeIdentifier getVolumeID() const {return m_volumeID;}
    /// Update cell identifier
    void setVolumeID( VolumeIdentifier value ) {m_volumeID = value;}
    /// Retrieve energy at this point
    Float_t getEnergy() const {return m_energy;}
    /// Set energy at this point
    void setEnergy(Float_t energy) {m_energy = energy;}
    /// Retrieve the hit point
    TVector3 getPoint() const {return m_point;}
    /// Set the hit point
    void setPoint(const TVector3& point) {m_point = point;}
    /// Overload the new operator to use our managed pool
    void* operator new(size_t size);
    /// Overload the new operator again for managed pool
    void* operator new(size_t size, void* vp);
    /// We need the = operator for pool management
    McTrajectoryPoint& operator=(const McTrajectoryPoint& rhs);

private:
    VolumeIdentifier m_volumeID;  // Volumen identifier of this point
    Float_t          m_energy;    // Trajectory total energy at this point
    TVector3         m_point;     // x,y,z coordinates of the point
    
    ClassDef(McTrajectoryPoint,1) // Monte Carlo Trajectory Class
};

/** @class McTrajectory
* @brief GLAST Monte Carlo Particle Trajectory class.
*
* The McTrajectory class contains the step by step position information for
* particles tracked through GLAST during the simulation stage. The class contains:
* - A reference to the McParticle to which it belongs
* - A TObjArray of positions as the particle traversed the GLAST universe
* - The charge of the particle 
* The MC data stored to ROOT files mirrors the data available within the Gaudi 
* TDS.
*  
*  $Header$
*/

class McTrajectory : public TObject  
{
public: 
    McTrajectory();
    
    virtual ~McTrajectory();
    
    /// clear lists, free pointers, etc., after read from / write to file
    void Clear(Option_t *option ="");
    void Fake( Int_t ievent, UInt_t rank, Float_t randNum ) ; // for tests
    Bool_t CompareInRange( const McTrajectory &, const std::string & name = "" ) const ; // for tests
    void Print(Option_t *option="") const;
    
    void initialize();

    /// Add McTrajectoryPoint to our list
    void addMcTrajectoryPoint(McTrajectoryPoint* point) {m_pointCol.Add(point);}

    /// return the full TObjArray containing the points on the trajectory
    const TObjArray* getMcPointCol() const { return &m_pointCol; };

    /// retrieve number of points stored in the collection
    inline UInt_t getMcPointCount() const 
    { 
//        return ((m_pointCol) ? m_pointCol->GetEntries() : 0);
        return m_pointCol.GetEntries();
    };
    
    /// add a TVector3 to points collection
    void addMcPoint(McTrajectoryPoint* point) {m_pointCol.Add(point);}

    /// clear of the array (necessary for converters)
    void clearMcPointCol()  { m_pointCol.Clear(); };

    /// Overload the new operator to use our managed pool
    void* operator new(size_t size);
    /// Overload the new operator again for managed pool
    void* operator new(size_t size, void* vp);
    /// We need the = operator for pool management
    McTrajectory& operator=(const McTrajectory& rhs);

private:
    /// An array of points defining the trajectory
    TObjArray m_pointCol;
    
    ClassDef(McTrajectory,2) // Monte Carlo Trajectory Class
};

#endif
