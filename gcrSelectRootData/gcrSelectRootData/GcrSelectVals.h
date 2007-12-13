#ifndef ROOT_GCRSELECTVALS_H
#define ROOT_GCRSELECTVALS_H

#include "TObject.h"
#include "TVector3.h"

/** @page gcrSelectRootData_doc
    @section GcrSelectVals GcrSelectVals Documentation
    Notes 

<table>
<tr><th><EM> Method <th> Return Type </EM><th> Description <th>

<small>
<tr><td> getInferedZ
<td>Int_t<td> Returns 
<tr><td> getAcdZ
<td>Int_t<td> Returns
<tr><td> getInteractionParmas
<td>Int_t<td> Returns 

</small>
</table>
*/

//-----------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------

class GcrSelectVals: public TObject
{
private:
    
    Int_t m_inferedZ;
    Int_t m_acdZ;
    Int_t m_interactionParams;

    
 public:
  
    GcrSelectVals()  {};
  
    ~GcrSelectVals()  {};
  
    GcrSelectVals(Int_t inferedZ, Int_t acdZ, Int_t interactionParams)
      :m_inferedZ(inferedZ),
       m_acdZ(acdZ),
       m_interactionParams(interactionParams)
       {}

    void Clear(Option_t *option="");
    void Print(Option_t *option="")const;

    void Fake(Int_t ievent, Float_t randNum);
    Bool_t CompareInRange(const GcrSelectVals &val, const std::string &name="") const;
  

    void                   setInferedZ (Int_t inferedZ)         {m_inferedZ  = inferedZ;}
    void                   setAcdZ  (Int_t acdZ)         {m_acdZ  = acdZ;}
    void                   setInteractionParams  (Int_t interactionParams)         {m_interactionParams  = interactionParams;}

    Int_t                    getInferedZ ()    const                      {return m_inferedZ    ;}
    Int_t                    getAcdZ ()        const                  {return m_acdZ    ;}
    Int_t                 getInteractionParams () const                          {return m_interactionParams    ;}

  
 
    ClassDef(GcrSelectVals,2)
};
    
#endif
  
