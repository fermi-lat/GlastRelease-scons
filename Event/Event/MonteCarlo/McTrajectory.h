#ifndef Event_McTrajectory_H
#define Event_McTrajectory_H 1

#include <iostream>
#include <vector>
#include <utility>
#include "CLHEP/Geometry/Point3D.h"
#include "GaudiKernel/Kernel.h"
#include "GaudiKernel/ContainedObject.h"
#include "GaudiKernel/SmartRefVector.h"
#include "Event/TopLevel/Definitions.h"
#include "idents/VolumeIdentifier.h"
#include "Event/Utilities/CLHEPStreams.h"
#include "Event/Utilities/IDStreams.h"
// Include all Glast container types here
//   to simplify inlude statements in algorithms
#include "GaudiKernel/ObjectVector.h"
#include "GaudiKernel/ObjectList.h"


/** @class McTrajectory
 * @brief Monte Carlo trajectory points
 *
 * @author:      R.Giannitrapani 
 *
 * $Header$
 */

#include "Event/MonteCarlo/McParticle.h"

extern const CLID& CLID_McTrajectory;

namespace Event {  // NameSpace

class McTrajectory : virtual public ContainedObject {
  public:

    virtual const CLID& clID() const   { return McTrajectory::classID(); }
    static const CLID& classID()       { return CLID_McTrajectory; }

    McTrajectory(){}
    ~McTrajectory(){}
    
    /// Add the 3d points to the trajectory
    void addPoints(std::vector<Hep3Vector>& points);
    /// Set the pointer to the McParticle
    void setMcParticle(SmartRef<McParticle> value);
    void setMcParticle( McParticle* value );
    /// Get the pointer to the McParticle
    const McParticle* getMcParticle() const;
    McParticle* getMcParticle();

    /// Get the 3d points
    std::vector<Hep3Vector>& getPoints(){return m_points;}

    /// Serialize the object for writing
    virtual StreamBuffer& serialize( StreamBuffer& s ) const;
    /// Serialize the object for reading
    virtual StreamBuffer& serialize( StreamBuffer& s );
    /// Fill the ASCII output stream
    virtual std::ostream& fillStream( std::ostream& s ) const;
  private:
    /// Pointer to McParticle of this trajectory
    SmartRef<McParticle>    m_mcParticle;
    /// The point of the trajectory
    std::vector<Hep3Vector> m_points;

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
