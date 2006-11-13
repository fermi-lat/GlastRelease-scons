
#include <RootConvert/Recon/GcrTrackConvert.h>
//#include <RootConvert/Recon/CalParamsConvert.h>
#include <RootConvert/Utilities/Toolkit.h>

namespace RootPersistence {

void convert( const Event::GcrTrack& gcrTrackTds, GcrTrack& gcrTrackRoot )
{

    //std::cout << "GcrTrackConvert::convert from TDS to ROOt BEGIN" << std::endl;   

	    Vector m_direction;
	    Vector m_dirError;
	    Point m_calEntryPoint;
	    Point m_calExitPoint;

     Vector directionTds = gcrTrackTds.getDirection();
     Vector dirErrorTds = gcrTrackTds.getDirError();
     Point calEntryPointTds = gcrTrackTds.getCalEntryPoint();
     Point calExitPointTds = gcrTrackTds.getCalExitPoint();


     TVector3 directionRoot(directionTds.x(), directionTds.y(), directionTds.z());
     TVector3 dirErrorRoot(dirErrorTds.x(), dirErrorTds.y(), dirErrorTds.z());
     TVector3 calEntryPointRoot(calEntryPointTds.x(), calEntryPointTds.y(), calEntryPointTds.z());
     TVector3 calExitPointRoot(calExitPointTds.x(), calExitPointTds.y(), calExitPointTds.z());
   
     gcrTrackRoot.setDirection(directionRoot);
     gcrTrackRoot.setDirError(dirErrorRoot);
     gcrTrackRoot.setCalEntryPoint(calEntryPointRoot);
     gcrTrackRoot.setCalExitPoint(calExitPointRoot);
        
   
    //  std::cout << "GcrTrackConvert::convert from TDS to ROOt END" << std::endl;   
    

}
 
void convert( const GcrTrack& gcrTrackRoot, Event::GcrTrack& gcrTrackTds )
{
    //std::cout << "GcrTrackConvert::convert from ROOT to TDS BEGIN" << std::endl;   
    
    
     TVector3 directionRoot=gcrTrackRoot.getDirection();
     TVector3 dirErrorRoot=gcrTrackRoot.getDirError();
     TVector3 calEntryPointRoot=gcrTrackRoot.getCalEntryPoint();
     TVector3 calExitPointRoot=gcrTrackRoot.getCalExitPoint();

    Vector   directionTds(directionRoot.x(),directionRoot.y(), directionRoot.z());
    Vector   dirErrorTds(dirErrorRoot.x(),dirErrorRoot.y(), dirErrorRoot.z());
    Point    calEntryPointTds(calEntryPointRoot.x(), calEntryPointRoot.y(), calEntryPointRoot.z());
    Point    calExitPointTds(calExitPointRoot.x(), calExitPointRoot.y(), calExitPointRoot.z());

    gcrTrackTds.setDirection(directionTds);
    gcrTrackTds.setDirError(dirErrorTds);
    gcrTrackTds.setCalEntryPoint(calEntryPointTds);
    gcrTrackTds.setCalExitPoint(calExitPointTds);

    
    //std::cout << "GcrTrackConvert::convert from ROOT to TDS END" << std::endl;   
 
  
}

}
