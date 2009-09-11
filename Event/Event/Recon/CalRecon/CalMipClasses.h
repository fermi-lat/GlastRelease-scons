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
namespace Event 
{
//-----------------------------------------------------------------------------------------------------------------
// Define a CalMipXtal class which will be used in the MIP finding
class CalMipXtal: virtual public ContainedObject 
{
private:
    Event::CalXtalRecData* m_xtalData;
    double                 m_d2C;
    bool                   m_free;
    bool                   m_freeC0;
    double                 m_ecor;

public:
    CalMipXtal() : m_free(true){};

    CalMipXtal(Event::CalXtalRecData* xtalData, double d2C, bool free, bool freeC0, double ecor) : 
                m_xtalData(xtalData), m_d2C(d2C), m_free(free), m_freeC0(freeC0), m_ecor(ecor) {};

    virtual ~CalMipXtal() {};

    void                   initialize(Event::CalXtalRecData* xtalData, double d2C, bool free, bool freC0, double ecor);

    void                   setD2C  (double d2C)                      {m_d2C      = d2C      ;}
    void                   setFree (bool free)                       {m_free     = free     ;}
    void                   setFreeC0 (bool freeC0)                   {m_freeC0   = freeC0   ;}
    void                   setXtal (Event::CalXtalRecData* xtalData) {m_xtalData = xtalData ;}
    void                   setEcor (double ene)                      {m_ecor     = ene;}

    double                 getD2C  ()                                {return m_d2C          ;}
    bool                   getFree ()                                {return m_free         ;}
    bool                   getFreeC0()                               {return m_freeC0       ;}
    Event::CalXtalRecData* getXtal ()                                {return m_xtalData     ;}
    double                 getEcor ()                                {return m_ecor         ;}

    /// Utilities 
    void writeOut(MsgStream& log) const; 
    std::ostream& fillStream( std::ostream& s ) const;
};

// Define a vector of CalMipXtals
typedef std::vector<CalMipXtal> CalMipXtalVec;

//-----------------------------------------------------------------------------------------------------------------
// Define a CalMipTrack class which will be used in the MIP finding
class CalMipTrack: public CalMipXtalVec, virtual public ContainedObject 
{
private:
    Point    m_point;
    Vector   m_dir;
    double   m_dt2C;
    double   m_d2Edge;
    int      m_calEdge;
    double   m_arcLen;
    double   m_ecor;
    double   m_ecorRms;
    double   m_chi2;
    double   m_erm;

public:

    CalMipTrack()  {};
  
    CalMipTrack(CalMipXtalVec calMipXtalTrack, Point point, Vector dir, double dt2C, double d2Edge, int calEdge, double arcLen, double ecor, double ecorRms, double chi2, double erm) :
      CalMipXtalVec(calMipXtalTrack), m_point(point), m_dir(dir), m_dt2C(dt2C), m_d2Edge(d2Edge), m_calEdge(calEdge), m_arcLen(arcLen), m_ecor(ecor), m_ecorRms(ecorRms), m_chi2(chi2), m_erm(erm)  {};
  
    ~CalMipTrack()  {};
  
    void initialize(CalMipXtalVec calMipXtalTrack, Point point, Vector dir, double dt2C, double d2Edge, int calEdge, double arcLen, double ecor, double ecorRms, double chi2, double erm);

    void    setPoint   (Point point)     {m_point    = point;    }
    void    setDir     (Vector dir)      {m_dir      = dir;      }
    void    setD2C     (double dt2C)     {m_dt2C     = dt2C;     }
    void    setD2Edge  (double d2Edge)   {m_d2Edge   = d2Edge;   }
    void    setCalEdge (int calEdge)     {m_calEdge  = calEdge;  }
    void    setArcLen  (double arcLen)   {m_arcLen   = arcLen;   }
    void    setEcor    (double ecor)     {m_ecor     = ecor;     }
    void    setEcorRms (double ecorRms)  {m_ecorRms  = ecorRms;  }
    void    setChi2    (double chi2)     {m_chi2     = chi2;     }
    void    setErm     (double erm)      {m_erm      = erm;      }

    int     getNh      ()    const       {return size();         }

    Point   getPoint   ()    const       {return m_point;        }
    Vector  getDir     ()    const       {return m_dir;          }
    double  getD2C     ()    const       {return m_dt2C;         }
    double  getD2Edge  ()    const       {return m_d2Edge;       }
    int     getCalEdge ()    const       {return m_calEdge;      }
    double  getArcLen  ()    const       {return m_arcLen;       }
    double  getEcor    ()    const       {return m_ecor;         }
    double  getEcorRms ()    const       {return m_ecorRms;      }
    double  getChi2    ()    const       {return m_chi2;         }
    double  getErm     ()    const       {return m_erm;          }

    /// Utilities 
    void writeOut(MsgStream& log) const; 
    std::ostream& fillStream( std::ostream& s ) const;
};

//-----------------------------------------------------------------------------------------------------------------
// Define a vector of CalMipTracks
typedef std::vector<CalMipTrack> CalMipTrackVec;

//typedef for the Container
typedef ObjectVector<CalMipTrack>      CalMipTrackCol;
typedef CalMipTrackCol::iterator       CalMipTrackColItr;
typedef CalMipTrackCol::const_iterator CalMipTrackColConItr;
}

#endif
