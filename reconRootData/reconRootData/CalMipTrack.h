#ifndef ROOT_CALMIPTRACK_H
#define ROOT_CALMIPTRACK_H

#include "TObject.h"
#include "TVector3.h"

#include "CalXtalRecData.h"

#include <vector>
#ifndef R__GLOBALSTL
#ifndef WIN32
using std::vector;
#else
using namespace std;
#endif
#endif

//-----------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------
// Define a CalMipTrack class which will be used in the MIP fitting
class CalMipTrack: public TObject
{
private:
    TVector3   m_point;
    TVector3   m_dir;
    Double_t   m_d2C;
    Double_t   m_d2Edge;
    Int_t      m_calEdge;
    Double_t   m_arcLen;
    Double_t   m_ecor;
    Double_t   m_ecorRms;
    Double_t   m_chi2;   
    Double_t   m_erm;   
  
 public:
  
    //CalMipXtalVec  m_calMipXtalTrack;
    //vector<CalXtalRecData> m_calMipXtalRecDataTrack;
  
    CalMipTrack()  {};
  
    ~CalMipTrack()  {};
  
    CalMipTrack(TVector3 point,
                TVector3 dir,
                Double_t d2C,
                Double_t d2Edge,
                Int_t    calEdge,
                Double_t arcLen,
                Double_t ecor,
                Double_t ecorRms,
                Double_t chi2,
                Double_t erm)
      :m_point(point),
       m_dir(dir),
       m_d2C(d2C),
       m_d2Edge(d2Edge),
       m_calEdge(calEdge),
       m_arcLen(arcLen),
       m_ecor(ecor),
       m_ecorRms(ecorRms),
       m_chi2(chi2),
       m_erm(erm){}
    
    void    setPoint   (TVector3 point)    {m_point    = point;    }
    void    setDir     (TVector3 dir)      {m_dir      = dir;      }
    void    setD2C     (Double_t d2C)      {m_d2C      = d2C;      }
    void    setD2Edge  (Double_t d2Edge)   {m_d2Edge   = d2Edge;   }
    void    setCalEdge (Int_t calEdge)     {m_calEdge  = calEdge;  }
    void    setArcLen  (Double_t arcLen)   {m_arcLen   = arcLen;   }
    void    setEcor    (Double_t ecor)     {m_ecor     = ecor;     }
    void    setEcorRms (Double_t ecorRms)  {m_ecorRms  = ecorRms;  }
    void    setChi2    (Double_t chi2)     {m_chi2     = chi2;     }
    void    setErm     (Double_t erm)      {m_erm      = erm;      }

    TVector3  getPoint   ()                {return m_point;        }
    TVector3  getDir     ()                {return m_dir;          }
    Double_t  getD2C     ()                {return m_d2C;          }
    Double_t  getD2Edge  ()                {return m_d2Edge;       }
    Int_t     getCalEdge ()                {return m_calEdge;      }
    Double_t  getArcLen  ()                {return m_arcLen;       }
    Double_t  getEcor    ()                {return m_ecor;         }
    Double_t  getEcorRms ()                {return m_ecorRms;      }
    Double_t  getChi2    ()                {return m_chi2;         }
    Double_t  getErm     ()                {return m_erm;          }
  
    void print();
  
    ClassDef(CalMipTrack,1)
};
    
#endif
  
