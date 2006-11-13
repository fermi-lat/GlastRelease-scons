#ifndef ROOT_GCRSELECTVALS_H
#define ROOT_GCRSELECTVALS_H

#include "TObject.h"
#include "TVector3.h"



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
  

    void                   setInferedZ (Int_t inferedZ)         {m_inferedZ  = inferedZ;}
    void                   setAcdZ  (Int_t acdZ)         {m_acdZ  = acdZ;}
    void                   setInteractionParams  (Int_t interactionParams)         {m_interactionParams  = interactionParams;}

    Int_t                    getInferedZ ()    const                      {return m_inferedZ    ;}
    Int_t                    getAcdZ ()        const                  {return m_acdZ    ;}
    Int_t                 getInteractionParams () const                          {return m_interactionParams    ;}

    //void print() const;
  
 
    ClassDef(GcrSelectVals,2)
};
    
#endif
  
