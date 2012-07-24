/** @file CalOverlay.h
*
* @author Tracy Usher
*
* $Header$
*/

#ifndef Event_CalOverlay_H
#define Event_CalOverlay_H 1

//include files
#include "GaudiKernel/ContainedObject.h"
#include "GaudiKernel/ObjectVector.h"
#include "GaudiKernel/IInterface.h"

#include "idents/CalXtalId.h"
#include "geometry/Point.h"

static const CLID& CLID_CalOverlay = InterfaceID("CalOverlay", 1, 0);

namespace Event {

/*!
* \class CalOverlay
* \brief TDS version of the Cal overlay data 
*
* Contains the information needed to overlay background data information on
* top of simulated data. 
*/
class CalOverlay : virtual public ContainedObject 
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
    CalOverlay() {};

    //! constructor with layer, tower and ToT, but with an empty list of strips
    CalOverlay(const idents::CalXtalId& calXtalId, const Point& p, const double& energy) : 
        m_calXtalId(calXtalId), m_position(p), m_energy(energy), m_status(0) 
    {
    };
    //! Destructor
    virtual ~CalOverlay() {
    };

    //! Retrieve pointer to class defininition structure
    virtual const CLID& clID() const   { return CalOverlay::classID(); }
    static const CLID& classID()       { return CLID_CalOverlay; }

    //! Retrieve xtal identifier
    idents::CalXtalId getCalXtalId() const {return m_calXtalId;}
    //! Retrieve the position
    const Point&      getPosition()  const {return m_position;}
    //! Retrive the energy
    double            getEnergy()    const {return m_energy;}

    //! Set/modify the crystal identifier
    void setCalXtalId(const idents::CalXtalId& xtalId) {m_calXtalId = idents::CalXtalId(xtalId);}
    //! Set the position
    void setPosition(const Point& position)            {m_position  = position;}
    //! Set the energy
    void setEnergy(const double& energy)               {m_energy    = energy;}

    //! Serialize the object for reading
    virtual StreamBuffer& serialize( StreamBuffer& s );
    //! Serialize the object for writing
    virtual StreamBuffer& serialize( StreamBuffer& s ) const;
    //! Fill the ASCII output stream
    virtual std::ostream& fillStream( std::ostream& s ) const;

    /// Retrieve the status word
    inline unsigned int getStatus() const {return m_status;}

    /// Allow status word to be modified
    inline void setStatus(unsigned int status)       {m_status  = status;}
    inline void addToStatus(DigiStatusBits bitToAdd) {m_status |= bitToAdd;}

private:

    idents::CalXtalId m_calXtalId;
    Point             m_position;
    double            m_energy;
    unsigned int      m_status;
};


//! Serialize the object for writing
inline StreamBuffer& CalOverlay::serialize( StreamBuffer& s ) const 
{
    ContainedObject::serialize(s);  
    s << m_calXtalId;

    return s;
}

//! Serialize the object for reading
inline StreamBuffer& CalOverlay::serialize( StreamBuffer& s )       
{
    ContainedObject::serialize(s);
    long calXtalId;
    s >> calXtalId;

    m_calXtalId = idents::CalXtalId(calXtalId);

    return s;
}

//! Fill the ASCII output stream

inline std::ostream& CalOverlay::fillStream( std::ostream& s ) const 
{
    s << "class CalOverlay :" << std::endl
      << "CalXtalId: " << m_calXtalId
      << std::endl;

    return s;
}

//! Definition of all container type of CalOverlay

typedef ObjectVector<CalOverlay> CalOverlayCol;

} //Namespace Event
#endif
