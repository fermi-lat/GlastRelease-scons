#include "gcrSelectRootData/GcrSelectVals.h"
#include <commonRootData/RootDataUtil.h>
#include "Riostream.h"

ClassImp(GcrSelectVals)


void GcrSelectVals::Clear(Option_t* /* option */)
{
    
    m_inferedZ = 0;
    m_acdZ = 0;
    m_interactionParams = 0;
}


void GcrSelectVals::Print(Option_t *option) const 
{
    TObject::Print(option);
    std::cout
        << "InferedZ: " << m_inferedZ << " acdZ: " << m_acdZ 
        << "interactionParams: " << m_interactionParams	
	<< std::endl;
}


//======================================================
// For Unit Tests
//======================================================

void GcrSelectVals::Fake( Int_t ievent, Float_t randNum ) {
    m_inferedZ = ievent * 1;
    m_acdZ = ievent *2;
    m_interactionParams = ievent * 3;
}

#define COMPARE_IN_RANGE(att) rootdatautil::CompareInRange(get ## att(),ref.get ## att(),#att)


Bool_t GcrSelectVals::CompareInRange(const GcrSelectVals & ref, const std::string & name )const  {

    bool result = true ;
    result = COMPARE_IN_RANGE(InferedZ) && result;
    result = COMPARE_IN_RANGE(AcdZ) && result;
    result = COMPARE_IN_RANGE(InteractionParams) && result;

    if (!result) {
        if ( name == "" ) {
            std::cout<<"Comparison ERROR for "<<ClassName()<<std::endl ;
        }
        else {
            std::cout<<"Comparison ERROR for "<<name<<std::endl ;
        }
    }
    return result ;

}
