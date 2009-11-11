#ifndef ROOT_GCRSELECTEDXTAL_H
#define ROOT_GCRSELECTEDXTAL_H

#include "TObject.h"
#include "TVector3.h"

#include "commonRootData/idents/CalXtalId.h"

/** @class GcrSelectedXtal
 * 
 * $Header$
*/

/** @page gcrSelectRootData_doc
    @section GcrSelectedXtal GcrSelectedXtal Documentation
    
    <P> A GcrSelectedXtal is a log crossed by the primary particle track, having an energy deposit > 100 MeV (if event passes HIP Filter) or > 5 MeV (if event passes Gamma, MIP or GCN filters), selected after rejection of nuclear interactions </P>

<table>
<tr>
<th>Method</th><th>Description</th><th>Return Value</th>
</tr>
<small>
<tr><td>  getXtalId ()</td><td>Returns the CalXtal Id</td><td>CalXtalId</td>
</tr>
<tr><td>  getRawEnergy ()</td><td>Returns the energy deposited in this Xtal</td><td>Double_t</td>
</tr>
<tr><td>  getPathLength ()</td><td>Returns the length of the track segment contained in this Xtal</td><td>Double_t</td>
</tr>
<tr><td>  getCorrEnergy ()</td><td>Returns the path-length corrected energy </td><td>Double_t</td>
</tr>
<tr><td>  getSelectGrade ()</td><td>Obsolete </td><td>Int_t</td>
</tr>
<tr><td>  getClosestFaceDist ()</td><td>Returns the distance to the closest Xtal face of the middle point of the track segment contained in that Xtal</td><td>Double_t</td>
</tr>
<tr><td>  getCrossedFaces ()</td><td>Returns a code that contains information about what Xtal faces were crossed by the particle</td><td>Int_t</td>
</tr>
<tr><td>  getEntryPoint ()</td><td>Returns the 3D coordinates of the entry point in the Xtal</td><td>TVector3</td>
</tr>
<tr><td>  getExitPoint ()</td><td>Returns the 3D coordinates of the exit point in the Xtal</td><td>TVector3</td>
</tr>


</small>
</table>

*/

//-----------------------------------------------------------------------------------------------------------------

class GcrSelectedXtal: public TObject
{
private:
    
    Int_t m_selectGrade;
    Int_t m_crossedFaces;
    CalXtalId m_xtalId;
    Double_t m_rawEnergy;
    Double_t m_pathLength; 
    Double_t m_corrEnergy;

    Double_t m_closestFaceDist;
    TVector3 m_entryPoint;
    TVector3 m_exitPoint;


    
 public:
  
    GcrSelectedXtal()  {};
  
    ~GcrSelectedXtal()  {};
  
    GcrSelectedXtal(CalXtalId xtalId, Double_t rawEnergy, Double_t pathLength, Double_t corrEnergy, Int_t selectGrade, Double_t closestFaceDist, Int_t crossedFaces, TVector3 entryPoint, TVector3 exitPoint)
      :m_selectGrade(selectGrade),
       m_crossedFaces(crossedFaces),
       m_xtalId(xtalId),
       m_rawEnergy(rawEnergy),
       m_pathLength(pathLength),
       m_corrEnergy(corrEnergy),
       m_closestFaceDist(closestFaceDist),
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



    void Print(Option_t *opt="") const;

    void print() const;

    void Fake(Int_t ievent, UInt_t ixtal, Float_t randNum);

    Bool_t CompareInRange( const GcrSelectedXtal & ref, const std::string & name="" ) const; 
  
 
    ClassDef(GcrSelectedXtal,5)
};
    
#endif
  
