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
    Double32_t m_d2C;
    Double32_t m_d2Edge;
    Int_t      m_calEdge;
    Double32_t m_arcLen;
    Double32_t m_ecor;
    Double32_t m_ecorRms;
    Double32_t m_chi2;   
    Double32_t m_erm;   
  
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

    TVector3  getPoint   () const          {return m_point;        }
    TVector3  getDir     () const          {return m_dir;          }
    Double_t  getD2C     () const          {return m_d2C;          }
    Double_t  getD2Edge  () const          {return m_d2Edge;       }
    Int_t     getCalEdge () const          {return m_calEdge;      }
    Double_t  getArcLen  () const          {return m_arcLen;       }
    Double_t  getEcor    () const          {return m_ecor;         }
    Double_t  getEcorRms () const          {return m_ecorRms;      }
    Double_t  getChi2    () const          {return m_chi2;         }
    Double_t  getErm     () const          {return m_erm;          }
  
    void print() const;
  
    void Fake( Int_t ievent, UInt_t rank, Float_t randNum ) ; // for tests
    Bool_t CompareInRange( const CalMipTrack &, const std::string & name = "" ) const ; // for tests

    ClassDef(CalMipTrack,3)
};
    
#endif
  
