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
    bool                   m_freeC0;

public:
    CalMipXtal() : m_free(true){};

    CalMipXtal(Event::CalXtalRecData* xtalData, double d2C, bool free, bool freeC0) : 
                m_xtalData(xtalData), m_d2C(d2C), m_free(free), m_freeC0(freeC0) {};

    ~CalMipXtal() {};

    void                   initialize(Event::CalXtalRecData* xtalData, double d2C, bool free, bool freeC0);

    void                   setD2C  (double d2C)                      {m_d2C      = d2C     ;}
    void                   setFree (bool free)                       {m_free     = free    ;}
    void                   setFreeC0 (bool freeC0)                       {m_freeC0     = freeC0    ;}
    void                   setXtal (Event::CalXtalRecData* xtalData) {m_xtalData = xtalData;}
    double                 getD2C  ()                                {return m_d2C         ;}
    bool                   getFree ()                                {return m_free        ;}
    bool                   getFreeC0()                                {return m_freeC0        ;}
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
    double   m_length;
    double   m_dTrack2C;
    double   m_dTrack2Edge;
    double   m_energy;

public:

    CalMipTrack()  {};
  
    CalMipTrack(Point point, Vector vector, int ndofTrack, double ki2Track, CalMipXtalVec calMipXtalTrack, double length, double dTrack2C, double dTrack2Edge, double energy) : 
             CalMipXtalVec(calMipXtalTrack), m_point(point), m_vector(vector), m_ndofTrack(ndofTrack), m_ki2Track(ki2Track),m_length(length), m_dTrack2C(dTrack2C), m_dTrack2Edge(dTrack2Edge), m_energy(energy)  {};
  
    ~CalMipTrack()  {};
  
    void initialize(Point point, Vector vector, int ndofTrack, double ki2Track, CalMipXtalVec calMipXtalTrack, double length, double dTrack2C, double dTrack2Edge, double energy);

    void           setPoint           (Point point)               {m_point = point;         }
    void           setDir             (Vector vector)             {m_vector = vector;       }
    void           setNdof            (int ndof)                  {m_ndofTrack = ndof;      }
    void           setKi2             (double ki2)                {m_ki2Track  = ki2;       }
    void           setLength          (double length)             {m_length      = length;     }
    void           setD2C             (double dTrack2C)           {m_dTrack2C    = dTrack2C;   }
    void           setD2Edge          (double dTrack2Edge)        {m_dTrack2Edge = dTrack2Edge;}
    void           setEnergy          (double energy)             {m_energy      = energy;     }

    Point          getPoint           ()    const                 {return m_point;          }
    Vector         getDir             ()    const                 {return m_vector;         }
    int            getNh              ()    const                 {return size();           }
    int            getNdof            ()    const                 {return m_ndofTrack;      }
    double         getKi2             ()    const                 {return m_ki2Track;       }
    double         getLength          ()    const                 {return m_length;            }
    double         getD2C             ()    const                 {return m_dTrack2C;          }
    double         getD2Edge          ()    const                 {return m_dTrack2Edge;       }
    double         getEnergy          ()    const                 {return m_energy;            }
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
