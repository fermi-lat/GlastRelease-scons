
#include <RootConvert/GcrSelect/GcrSelectValsConvert.h>
//#include <RootConvert/Recon/CalParamsConvert.h>
#include <RootConvert/Utilities/Toolkit.h>
#include <cmath>

namespace RootPersistence {

void convert( const Event::GcrSelectVals& gcrSelectValsTds, GcrSelectVals& gcrSelectValsRoot )
{

    //std::cout << "GcrSelectValsConvert::convert from TDS to ROOt BEGIN" << std::endl; 
    
    
    Int_t inferedZ=gcrSelectValsTds.getInferedZ();
    Int_t acdZ=gcrSelectValsTds.getAcdZ();
    Int_t interactionParams=(Int_t)floor(gcrSelectValsTds.getInteractionParams());
    UInt_t gcrOBFStatusWord=gcrSelectValsTds.getGcrOBFStatusWord();

    gcrSelectValsRoot.setInferedZ(inferedZ);
    gcrSelectValsRoot.setAcdZ(acdZ);
    gcrSelectValsRoot.setInteractionParams(interactionParams);
    gcrSelectValsRoot.setGcrOBFStatusWord(gcrOBFStatusWord);
   
 
   
    //  std::cout << "GcrSelectValsConvert::convert from TDS to ROOt END" << std::endl;   
    

}
 
void convert( const GcrSelectVals& gcrSelectValsRoot, Event::GcrSelectVals& gcrSelectValsTds )
{
    //std::cout << "GcrSelectValsConvert::convert from ROOT to TDS BEGIN" << std::endl;   
    
    
    Int_t inferedZ=gcrSelectValsRoot.getInferedZ();
    Int_t acdZ=gcrSelectValsRoot.getAcdZ();
    Int_t interactionParams=gcrSelectValsRoot.getInteractionParams();
    UInt_t gcrOBFStatusWord=gcrSelectValsRoot.getGcrOBFStatusWord();

    gcrSelectValsTds.setInferedZ(inferedZ);
    gcrSelectValsTds.setAcdZ(acdZ);
    gcrSelectValsTds.setInteractionParams(interactionParams);
    gcrSelectValsTds.setGcrOBFStatusWord(gcrOBFStatusWord);
 
 
    
    //std::cout << "GcrSelectValsConvert::convert from ROOT to TDS END" << std::endl;   
 
  
}

}
