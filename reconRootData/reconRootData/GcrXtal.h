#ifndef ROOT_GCRXTAL_H
#define ROOT_GCRXTAL_H

#include "TObject.h"
#include "TVector3.h"

//#include "CalXtalRecData.h"
#include "commonRootData/idents/CalXtalId.h"


class GcrXtal: public TObject
{
private:
    
    CalXtalId m_xtalId;
    Double32_t m_pathLength; 
    Double32_t m_closestFaceDist;
    Int_t m_crossedFaces;
    TVector3 m_entryPoint;
    TVector3 m_exitPoint;
    
    
 public:
  
    GcrXtal()  {};
  
    ~GcrXtal()  {};
  
    GcrXtal(CalXtalId xtalId, Double_t pathLength, Double_t closestFaceDist, Int_t crossedFaces, TVector3 entryPoint, TVector3 exitPoint)
      :m_xtalId(xtalId),
       m_pathLength(pathLength),
       m_closestFaceDist(closestFaceDist),
       m_crossedFaces(crossedFaces),
       m_entryPoint(entryPoint),
       m_exitPoint(exitPoint)
       {}
  
    void    setXtalId (CalXtalId xtalId) {m_xtalId = xtalId ;}
    void    setPathLength  (Double_t pathLength)         {m_pathLength  = pathLength;}

    void                   setClosestFaceDist (Double_t closestFaceDist) {m_closestFaceDist = closestFaceDist ;}
    void                   setCrossedFaces  (Int_t crossedFaces)         {m_crossedFaces  = crossedFaces;}
    void                   setEntryPoint (TVector3 entryPoint) {m_entryPoint = entryPoint ;}
    void                   setExitPoint  (TVector3 exitPoint)         {m_exitPoint  = exitPoint;}


    CalXtalId getXtalId ()       const                          {return m_xtalId     ;}
    Double_t         getPathLength () const                         {return m_pathLength    ;}

    Double_t                  getClosestFaceDist () const                          {return m_closestFaceDist    ;}
    Int_t                  getCrossedFaces () const                          {return m_crossedFaces    ;}
    TVector3                  getEntryPoint () const                          {return m_entryPoint    ;}
    TVector3                  getExitPoint () const                          {return m_exitPoint    ;}


    void Print(Option_t *option="") const;

 
    void Fake( Int_t ievent, Int_t ixtal, Float_t randNum ) ; // for tests
    Bool_t CompareInRange(const GcrXtal&, const std::string & name = "" )const ; // for tests


    void print() const;
  
 
    ClassDef(GcrXtal,5)
};
    
#endif
  
