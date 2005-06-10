#ifndef CalMipClasses_H
#define CalMipClasses_H

#include "geometry/Point.h"
#include "geometry/Vector.h"

#include "Event/Recon/CalRecon/CalXtalRecData.h"

#include "GaudiKernel/ObjectVector.h"
#include "GaudiKernel/ContainedObject.h"
#include "GaudiKernel/MsgStream.h"

static const CLID& CLID_CalMipTrackVecCol = InterfaceID("CalMipTrackVecCol", 1, 0);

/**   
* @class StdMipFindingTool
*
* $Header$
*/


//-----------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------
// Define a CalMipXtal class which will be used in the MIP fitting

namespace Event 
{
//-----------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------
class CalMipXtal: virtual public ContainedObject 
{
private:
    Event::CalXtalRecData* m_xtalData;
    double                 m_d2C;
    bool                   m_free;

public:
    CalMipXtal() : m_free(true){};

    CalMipXtal(Event::CalXtalRecData* xtalData, double d2C, bool free) : 
                m_xtalData(xtalData), m_d2C(d2C), m_free(free) {};

    ~CalMipXtal() {};

    void                   initialize(Event::CalXtalRecData* xtalData, double d2C, bool free);

    void                   setD2C  (double d2C)                      {m_d2C      = d2C     ;}
    void                   setFree (bool free)                       {m_free     = free    ;}
    void                   setXtal (Event::CalXtalRecData* xtalData) {m_xtalData = xtalData;}
    double                 getD2C  ()                                {return m_d2C         ;}
    bool                   getFree ()                                {return m_free        ;}
    Event::CalXtalRecData* getXtal ()                                {return m_xtalData    ;}

    /// Utilities 
    void writeOut(MsgStream& log) const; 
    std::ostream& fillStream( std::ostream& s ) const;
};

// Define a vector of CalMipXtal
typedef std::vector<CalMipXtal> CalMipXtalVec;


//-----------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------
// Define a CalMipTrack class which will be used in the MIP fitting
class CalMipTrack: public CalMipXtalVec, virtual public ContainedObject 
{
private:
    Point    m_point;
    Vector   m_vector;
    int      m_ndofTrack;
    double   m_ki2Track;

public:

    CalMipTrack()  {};
  
    CalMipTrack(Point point, Vector vector, int ndofTrack, double ki2Track, CalMipXtalVec calMipXtalTrack) : 
             CalMipXtalVec(calMipXtalTrack), m_point(point), m_vector(vector), m_ndofTrack(ndofTrack), m_ki2Track(ki2Track) {};
  
    ~CalMipTrack()  {};
  
    void initialize(Point point, Vector vector, int ndofTrack, double ki2Track, CalMipXtalVec calMipXtalTrack);

    void           setPoint           (Point point)               {m_point = point;         }
    void           setDir             (Vector vector)             {m_vector = vector;       }
    void           setNdof            (int ndof)                  {m_ndofTrack = ndof;      }
    void           setKi2             (double ki2)                {m_ki2Track  = ki2;       }

    Point          getPoint           ()    const                 {return m_point;          }
    Vector         getDir             ()    const                 {return m_vector;         }
    int            getNh              ()    const                 {return size();           }
    int            getNdof            ()    const                 {return m_ndofTrack;      }
    double         getKi2             ()    const                 {return m_ki2Track;       }

    /// Utilities 
    void writeOut(MsgStream& log) const; 
    std::ostream& fillStream( std::ostream& s ) const;
};

//-----------------------------------------------------------------------------------------------------------------
// Define a vector of CalMipXtal
typedef std::vector<CalMipTrack> CalMipTrackVec;

//typedef for the Container
typedef ObjectVector<CalMipTrack>      CalMipTrackCol;
typedef CalMipTrackCol::iterator       CalMipTrackColItr;
typedef CalMipTrackCol::const_iterator CalMipTrackColConItr;
}

#endif
