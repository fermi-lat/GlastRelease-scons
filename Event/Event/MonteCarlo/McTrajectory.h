#ifndef Event_McTrajectory_H
#define Event_McTrajectory_H 1

#include <iostream>
#include <vector>
#include <map>
#include <utility>
#include "CLHEP/Geometry/Point3D.h"
#include "GaudiKernel/Kernel.h"
#include "GaudiKernel/ContainedObject.h"
#include "GaudiKernel/SmartRefVector.h"
#include "Event/TopLevel/Definitions.h"
#include "idents/VolumeIdentifier.h"
#include "Event/Utilities/CLHEPStreams.h"
#include "Event/Utilities/IDStreams.h"
#include "Event/MonteCarlo/McParticle.h"
#include "Event/MonteCarlo/McPositionHit.h"
#include "Event/MonteCarlo/McIntegratingHit.h"
// Include all Glast container types here
//   to simplify inlude statements in algorithms
#include "GaudiKernel/ObjectVector.h"
#include "GaudiKernel/ObjectList.h"

#include "idents/VolumeIdentifier.h"

/** @class McTrajectory
 * @brief Monte Carlo trajectory points
 *
 * @author:      R.Giannitrapani 
 *
 * $Header$
 */

#include "GaudiKernel/IInterface.h"

static const CLID& CLID_McTrajectory = InterfaceID("McTrajectory", 1, 0);

namespace Event {  // NameSpace

// Class to hold information at each point on the trajectory
class McTrajectoryPoint
{
public:
    McTrajectoryPoint() : m_energy(0.), m_point(0.,0.,0.) {}
    McTrajectoryPoint(idents::VolumeIdentifier vId, float energy, CLHEP::Hep3Vector& point) :
        m_volumeID(vId), m_energy(energy), m_point(point) {}
    virtual ~McTrajectoryPoint() {}

    /// Retrieve cell identifier
    idents::VolumeIdentifier getVolumeID() const {return m_volumeID;}
    /// Update cell identifier
    void setVolumeID( idents::VolumeIdentifier value ) {m_volumeID = value;}
    /// Retrieve energy at this point
    float getEnergy() const {return m_energy;}
    /// Set energy at this point
    void setEnergy(float energy) {m_energy = energy;}
    /// Retrieve the hit point
    CLHEP::Hep3Vector getPoint() const {return m_point;}
    /// Set the hit point
    void setPoint(const CLHEP::Hep3Vector& point) {m_point = point;}
    /// Fill the ASCII output stream
    virtual std::ostream& fillStream( std::ostream& s ) const {return s << "energy=" << m_energy << " ";}
        
    friend std::ostream& operator << (std::ostream& s, const McTrajectoryPoint& obj)
    {
        return obj.fillStream(s);
    };

private:
    idents::VolumeIdentifier m_volumeID;  // Volumen identifier of this point
    float                    m_energy;    // Trajectory total energy at this point
    CLHEP::Hep3Vector        m_point;     // x,y,z coordinates of the point
};

class McTrajectory : virtual public ContainedObject 
{
public:

    virtual const CLID& clID() const   { return McTrajectory::classID(); }
    static const CLID& classID()       { return CLID_McTrajectory; }

    McTrajectory(){}
    ~McTrajectory();
    
    /// Add the 3d points to the trajectory
    void addPoints(std::vector<Event::McTrajectoryPoint*>& points);
    /// Add a 3d point to the collection
    void addPoint(Event::McTrajectoryPoint* point);
    /// Set the pointer to the McParticle
    //void setMcParticle(SmartRef<McParticle> value);
    //void setMcParticle( McParticle* value );
    /// Get the pointer to the McParticle
    //const McParticle* getMcParticle() const;
    //McParticle* getMcParticle();

    /// Get the 3d points
    const std::vector<Event::McTrajectoryPoint*>& getPoints() const {return m_points;}

    /// get, set charge
    //int getCharge() const { return m_charge; }
    //void setCharge(int charge){ m_charge=charge;}

    /// Serialize the object for writing
    virtual StreamBuffer& serialize( StreamBuffer& s ) const;
    /// Serialize the object for reading
    virtual StreamBuffer& serialize( StreamBuffer& s );
    /// Fill the ASCII output stream
    virtual std::ostream& fillStream( std::ostream& s ) const;
        
    friend std::ostream& operator << (std::ostream& s, const McTrajectory& obj)
    {
        return obj.fillStream(s);
    };

private:
    /// Pointer to McParticle of this trajectory
    //SmartRef<McParticle>    m_mcParticle;
    /// The point of the trajectory
    std::vector<McTrajectoryPoint*> m_points;
    /// the (redundant?) charge
    //int m_charge;

};

/// Serialize the object for writing
inline StreamBuffer& McTrajectory::serialize( StreamBuffer& s ) const
  {
    ContainedObject::serialize(s);
    return s
      << m_points.size();
  }


/// Serialize the object for reading
inline StreamBuffer& McTrajectory::serialize( StreamBuffer& s )
{
    ContainedObject::serialize(s);
    return s;
}


  

/// Fill the ASCII output stream
inline std::ostream& McTrajectory::fillStream( std::ostream& s ) const
{
    s << "class McTrajectory :"
      << "\n    Total number of points        = "
      << m_points.size();
    return s;
}


// Definition of all container types of McTrajectory
//template <class TYPE> class ObjectVector;
typedef ObjectVector<McTrajectory>     McTrajectoryVector;
typedef ObjectVector<McTrajectory>     McTrajectoryCol;
//template <class TYPE> class ObjectList;
typedef ObjectList<McTrajectory>       McTrajectoryList;


}

#endif // Event_McTrajectory_H
