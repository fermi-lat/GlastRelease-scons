
#include <RootConvert/GcrSelect/GcrSelectValsConvert.h>
//#include <RootConvert/Recon/CalParamsConvert.h>
#include <RootConvert/Utilities/Toolkit.h>

namespace RootPersistence {

void convert( const Event::GcrSelectVals& gcrSelectValsTds, GcrSelectVals& gcrSelectValsRoot )
{

    //std::cout << "GcrSelectValsConvert::convert from TDS to ROOt BEGIN" << std::endl; 
    
    
    Int_t inferedZ=gcrSelectValsTds.getInferedZ();
    Int_t acdZ=gcrSelectValsTds.getAcdZ();
    Int_t interactionParams=gcrSelectValsTds.getInteractionParams();

    gcrSelectValsRoot.setInferedZ(inferedZ);
    gcrSelectValsRoot.setAcdZ(acdZ);
    gcrSelectValsRoot.setInteractionParams(interactionParams);
   
 
   
    //  std::cout << "GcrSelectValsConvert::convert from TDS to ROOt END" << std::endl;   
    

}
 
void convert( const GcrSelectVals& gcrSelectValsRoot, Event::GcrSelectVals& gcrSelectValsTds )
{
    //std::cout << "GcrSelectValsConvert::convert from ROOT to TDS BEGIN" << std::endl;   
    
    
    Int_t inferedZ=gcrSelectValsRoot.getInferedZ();
    Int_t acdZ=gcrSelectValsRoot.getAcdZ();
    Int_t interactionParams=gcrSelectValsRoot.getInteractionParams();

    gcrSelectValsTds.setInferedZ(inferedZ);
    gcrSelectValsTds.setAcdZ(acdZ);
    gcrSelectValsTds.setInteractionParams(interactionParams);
 
    
    //std::cout << "GcrSelectValsConvert::convert from ROOT to TDS END" << std::endl;   
 
  
}

}
