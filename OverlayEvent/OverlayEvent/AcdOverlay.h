/** @file AcdOverlay.h
*
* @author Tracy Usher
*
* $Header$
*/

#ifndef Event_AcdOverlay_H
#define Event_AcdOverlay_H 1

//include files
#include "GaudiKernel/ContainedObject.h"
#include "GaudiKernel/ObjectVector.h"

#include "CLHEP/Geometry/Point3D.h"

#include "idents/AcdId.h"

static const CLID& CLID_AcdOverlay = InterfaceID("AcdOverlay", 1, 0);

namespace Event {

/*!
* \class AcdOverlay
* \brief TDS version of the Cal overlay data 
*
* Contains the information needed to overlay background data information on
* top of simulated data. 
*/
class AcdOverlay : virtual public ContainedObject 
{
public:
    //! typedefs
    typedef std::vector<int> HitList;
    typedef HitList::const_iterator const_iterator;

    /// Status bit definitions
    enum DigiStatusBits {
         DIGI_NULL    = 0,
         DIGI_MC      = 1<<28,   // This digi is from Monte Carlo
         DIGI_DATA    = 1<<30,   // This digi is from data 
         DIGI_OVERLAY = 1<<31    // This digi contains overlay info
    };

    //! Constructors
    //! Null constructor
    AcdOverlay() {};

    //! constructor with layer, tower and ToT, but with an empty list of strips
    AcdOverlay(const idents::VolumeIdentifier& volId, const idents::AcdId& acdId,
               const double energyDep, const HepPoint3D& point) : 
        m_volumeId(volId), m_acdId(acdId), m_energyDep(energyDep), m_center(point), m_status(0) 
    {
    };
    //! Destructor
    virtual ~AcdOverlay() {
    };

    //! Retrieve pointer to class defininition structure
    virtual const CLID& clID() const   { return AcdOverlay::classID(); }
    static const CLID& classID()       { return CLID_AcdOverlay; }

    //! Retrive the volume identifier
    const idents::VolumeIdentifier& getVolumeId()  const {return m_volumeId;}
    //! Retrieve xtal identifier
    const idents::AcdId&            getAcdId()     const {return m_acdId;}
    //! Retrieve the deposited energy
    const double                    getEnergyDep() const {return m_energyDep;}
    //! Retrieve the coordinates 
    const HepPoint3D&               getPosition()  const {return m_center;}
    //! Retrieve the status word
    inline unsigned int             getStatus()    const {return m_status;}

    //! Serialize the object for reading
    virtual StreamBuffer& serialize( StreamBuffer& s );
    //! Serialize the object for writing
    virtual StreamBuffer& serialize( StreamBuffer& s ) const;
    //! Fill the ASCII output stream
    virtual std::ostream& fillStream( std::ostream& s ) const;

    /// Allow status word to be modified
    inline void setStatus(unsigned int status)       {m_status  = status;}
    inline void addToStatus(DigiStatusBits bitToAdd) {m_status |= bitToAdd;}

    /// Set routines for the rest of the data members
    inline void setVolumeId(const idents::VolumeIdentifier& volId) {m_volumeId = volId;}
    inline void setAcdId(const idents::AcdId& acdId)               {m_acdId = acdId;}
    inline void setEnergyDep(double energy)                        {m_energyDep = energy;}
    inline void setPosition(const HepPoint3D& point)               {m_center = point;}

private:

    idents::VolumeIdentifier m_volumeId;
    idents::AcdId            m_acdId;
    double                   m_energyDep;
    HepPoint3D               m_center;
    unsigned int             m_status;
};


//! Serialize the object for writing
inline StreamBuffer& AcdOverlay::serialize( StreamBuffer& s ) const 
{
    ContainedObject::serialize(s);  
    s << m_acdId;

    return s;
}

//! Serialize the object for reading
inline StreamBuffer& AcdOverlay::serialize( StreamBuffer& s )       
{
    ContainedObject::serialize(s);
    long AcdId;
    s >> AcdId;

    m_acdId = idents::AcdId(AcdId);

    return s;
}

//! Fill the ASCII output stream

inline std::ostream& AcdOverlay::fillStream( std::ostream& s ) const 
{
    s << "class AcdOverlay :" << std::endl
      << "AcdId: " << m_acdId
      << std::endl;

    return s;
}

//! Definition of all container type of AcdOverlay

typedef ObjectVector<AcdOverlay> AcdOverlayCol;

} //Namespace Event
#endif
