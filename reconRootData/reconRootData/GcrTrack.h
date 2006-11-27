#ifndef ROOT_GCRTRACK_H
#define ROOT_GCRTRACK_H

#include "TObject.h"
#include "TVector3.h"

#include "CalXtalRecData.h"
#include "commonRootData/idents/CalXtalId.h"


//#include <vector>
//#include <point>
#ifndef R__GLOBALSTL
#ifndef WIN32
//using std::vector;
//using std::point;
#else
using namespace std;

#endif
#endif

//-----------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------

class GcrTrack: public TObject
{
private:
    
    TVector3 m_direction;
    TVector3 m_dirError;
    TVector3 m_calEntryPoint;
    TVector3 m_calExitPoint;

    
 public:
  
    GcrTrack()  {};
  
    ~GcrTrack()  {};
  
    GcrTrack(TVector3 direction, TVector3 dirError, TVector3 calEntryPoint, TVector3 calExitPoint)
      :m_direction(direction),
       m_dirError(dirError),
       m_calEntryPoint(calEntryPoint),
       m_calExitPoint(calExitPoint)
       {}
  

    void                   setDirection (TVector3 direction)         {m_direction  = direction;}
    void                   setDirError  (TVector3 dirError)         {m_dirError  = dirError;}
    void                   setCalEntryPoint  (TVector3 calEntryPoint)         {m_calEntryPoint  = calEntryPoint;}
    void                   setCalExitPoint  (TVector3 calExitPoint)         {m_calExitPoint  = calExitPoint;}

    TVector3                  getDirection ()    const                      {return m_direction    ;}
    TVector3                  getDirError ()     const                     {return m_dirError    ;}
    TVector3                  getCalEntryPoint ()   const                       {return m_calEntryPoint    ;}
    TVector3                  getCalExitPoint ()    const                      {return m_calExitPoint    ;}

    void Print(Option_t *option="") const;
  
 
    ClassDef(GcrTrack,2)
};
    
#endif
  
