#ifndef ROOT_GCRSELECTEDXTAL_H
#define ROOT_GCRSELECTEDXTAL_H

#include "TObject.h"
#include "TVector3.h"

#include <commonRootData/idents/CalXtalId.h>


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

class GcrSelectedXtal: public TObject
{
private:
    
    CalXtalId m_xtalId;
    Double_t m_rawEnergy;
    Double_t m_pathLength; 
    Double_t m_corrEnergy;
    Int_t m_selectGrade;

    Double_t m_closestFaceDist;
    Int_t m_crossedFaces;
    TVector3 m_entryPoint;
    TVector3 m_exitPoint;


    
 public:
  
    GcrSelectedXtal()  {};
  
    ~GcrSelectedXtal()  {};
  
    GcrSelectedXtal(CalXtalId xtalId, Double_t rawEnergy, Double_t pathLength, Double_t corrEnergy, Int_t selectGrade, Double_t closestFaceDist, Int_t crossedFaces, TVector3 entryPoint, TVector3 exitPoint)
      :m_xtalId(xtalId),
       m_rawEnergy(rawEnergy),
       m_pathLength(pathLength),
       m_corrEnergy(corrEnergy),
       m_selectGrade(selectGrade),
       m_closestFaceDist(closestFaceDist),
       m_crossedFaces(crossedFaces),
       m_entryPoint(entryPoint),
       m_exitPoint(exitPoint)
       {}
  
    void    setXtalId (CalXtalId xtalId) {m_xtalId = xtalId ;}
    void    setRawEnergy  (Double_t rawEnergy)         {m_rawEnergy  = rawEnergy;}
    void    setPathLength  (Double_t pathLength)         {m_pathLength  = pathLength;}
    void    setCorrEnergy  (Double_t corrEnergy)         {m_corrEnergy  = corrEnergy;}
    void    setSelectGrade  (Int_t selectGrade)         {m_selectGrade  = selectGrade;}
 
    void                   setClosestFaceDist (Double_t closestFaceDist) {m_closestFaceDist = closestFaceDist ;}
    void                   setCrossedFaces  (Int_t crossedFaces)         {m_crossedFaces  = crossedFaces;}
    void                   setEntryPoint (TVector3 entryPoint) {m_entryPoint = entryPoint ;}
    void                   setExitPoint  (TVector3 exitPoint)         {m_exitPoint  = exitPoint;}

    CalXtalId getXtalId ()              const                  {return m_xtalId     ;}
    Double_t         getRawEnergy ()    const                      {return m_rawEnergy    ;}
    Double_t         getPathLength ()   const                       {return m_pathLength    ;}
    Double_t         getCorrEnergy ()   const                       {return m_corrEnergy    ;}
    Int_t           getSelectGrade ()   const                       {return m_selectGrade    ;}

    Double_t                  getClosestFaceDist ()    const                      {return m_closestFaceDist    ;}
    Int_t                  getCrossedFaces ()         const                 {return m_crossedFaces    ;}
    TVector3                  getEntryPoint ()        const                  {return m_entryPoint    ;}
    TVector3                  getExitPoint ()         const                 {return m_exitPoint    ;}

    /// suggested by Zach: methods allowing to break variable crossedFaces comming soon... 
    typedef enum {
        XFACE_ZTOP,
        XFACE_ZBOT,
        XFACE_XLEFT,
        XFACE_XRIGHT,
        XFACE_YLEFT,
        XFACE_YRIGHT
    } XFACE_BITPOS;



    void print() const;
  
 
    ClassDef(GcrSelectedXtal,4)
};
    
#endif
  
