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
    TVector3          m_point;
    TVector3          m_dir;
    Int_t             m_ndofTrack;
    Double_t          m_ki2Track;
    Double_t          m_length;
    Double_t          m_dTrack2C;
    Double_t          m_dTrack2Edge;
    Double_t          m_energy;
  
 public:
  
    //CalMipXtalVec  m_calMipXtalTrack;
    //vector<CalXtalRecData> m_calMipXtalRecDataTrack;
  
    CalMipTrack()  {};
  
    ~CalMipTrack()  {};
  
    CalMipTrack(TVector3 point, 
                TVector3 dir, 
                Int_t    ndofTrack, 
                Double_t ki2Track, 
                Double_t length, 
                Double_t dTrack2C, 
                Double_t dTrack2Edge, 
                Double_t energy) 
                : m_point(point),
                  m_dir(dir),
                  m_ndofTrack(ndofTrack),
                  m_ki2Track(ki2Track),
                  m_length(length),
                  m_dTrack2C(dTrack2C),
                  m_dTrack2Edge(dTrack2Edge),
                  m_energy(energy){}
  
    void      setPoint(TVector3 point)        {m_point       = point;      }
    void      setDir(TVector3 vector)         {m_dir         = vector;     }
    void      setNdof(Int_t ndof)             {m_ndofTrack   = ndof;       }
    void      setKi2(Double_t ki2)            {m_ki2Track    = ki2;        }
    void      setLength(Double_t length)      {m_length      = length;     }
    void      setD2C(Double_t dTrack2C)       {m_dTrack2C    = dTrack2C;   }
    void      setD2Edge(Double_t dTrack2Edge) {m_dTrack2Edge = dTrack2Edge;}
    void      setEnergy(Double_t energy)      {m_energy      = energy;     }
    TVector3  getPoint()                      {return m_point;             }
    TVector3  getDir()                        {return m_dir;               }
    Int_t     getNdof()                       {return m_ndofTrack;         }
    Double_t  getKi2()                        {return m_ki2Track;          }
    Double_t  getLength()                     {return m_length;            }
    Double_t  getD2C()                        {return m_dTrack2C;          }
    Double_t  getD2Edge()                     {return m_dTrack2Edge;       }
    Double_t  getEnergy()                     {return m_energy;            }
  
    void print();
  
    ClassDef(CalMipTrack,1)
};
    
#endif
  
