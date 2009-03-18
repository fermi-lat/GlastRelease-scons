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
         PMTA_ACCEPT_BIT              = 0x1,    // channel is above zero suppresion threshold (applied to digital data)
         PMTA_VETO_BIT                = 0x2,    // channel fired veto discriminator (applied to analog data)
         PMTA_RANGE_BIT               = 0x4,    // channel was read out in high range
         PMTA_CNO_BIT                 = 0x8,    // could be any channel on that GARC
         PMTA_ODD_PARITY_ERROR_BIT    = 0x10,   // PHA data transmission had parity error
         PMTA_HEADER_PARITY_ERROR_BIT = 0x20,   // header data transmission had parity error 
         PMTA_DEAD_BIT                = 0x40,   // PMT was flagged as dead by offline calibration
         PMTA_HOT_BIT                 = 0x80,   // PMT was flagged as hot by offline calibration
         PMTA_ANY_ERROR_MASK          = 0xF0,   // PMT has any error
         PMTB_ACCEPT_BIT              = 0x100,  // channel is above zero suppresion threshold (applied to digital data)
         PMTB_VETO_BIT                = 0x200,  // channel fired veto discriminator (applied to analog data)
         PMTB_RANGE_BIT               = 0x400,  // channel was read out in high range
         PMTB_CNO_BIT                 = 0x800,  // could be any channel on that GARC
         PMTB_ODD_PARITY_ERROR_BIT    = 0x1000, // PHA data transmission had parity error
         PMTB_HEADER_PARITY_ERROR_BIT = 0x2000, // header data transmission had parity error 
         PMTB_DEAD_BIT                = 0x4000, // PMT was flagged as dead by offline calibration
         PMTB_HOT_BIT                 = 0x8000, // PMT was flagged as hot by offline calibration
         PMTB_ANY_ERROR_MASK          = 0xF000  // PMT has any error
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
