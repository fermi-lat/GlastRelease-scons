
#include <RootConvert/Recon/GcrXtalConvert.h>
#include <RootConvert/Recon/CalXtalRecDataConvert.h>
//#include <RootConvert/Recon/CalParamsConvert.h>
#include <RootConvert/Utilities/Toolkit.h>

namespace RootPersistence {

void convert( const Event::GcrXtal& gcrXtalTds, GcrXtal& gcrXtalRoot )
{

    //std::cout << "GcrXtalConvert::convert from TDS to ROOt BEGIN" << std::endl;   


    
    //Event::CalXtalRecData* xTalTds = gcrXtalTds.getXtal ();
    
    CalXtalId rootId = RootPersistence::convert(gcrXtalTds.getXtalId()) ; 
    double pathLength =   gcrXtalTds.getPathLength();
    double closestFaceDist = gcrXtalTds.getClosestFaceDist();
    int crossedFaces = gcrXtalTds.getCrossedFaces();  
   
    
   // CalXtalRecData* xTalRoot = new CalXtalRecData();
   // convert(*xTalTds, *xTalRoot);
    
    gcrXtalRoot.setXtalId(rootId);
    gcrXtalRoot.setPathLength(pathLength);
    gcrXtalRoot.setClosestFaceDist(closestFaceDist);
    gcrXtalRoot.setCrossedFaces(crossedFaces);
    
    Point entryPointTds = gcrXtalTds.getEntryPoint();
    Point exitPointTds = gcrXtalTds.getExitPoint();
    
    TVector3 entryPointRoot(entryPointTds.x(), entryPointTds.y(), entryPointTds.z());
    TVector3 exitPointRoot(exitPointTds.x(), exitPointTds.y(), exitPointTds.z());
    
    gcrXtalRoot.setEntryPoint(entryPointRoot);
    gcrXtalRoot.setExitPoint(exitPointRoot);
    
    //  std::cout << "GcrXtalConvert::convert from TDS to ROOt END" << std::endl;   
    

}
 
void convert( const GcrXtal& gcrXtalRoot, Event::GcrXtal& gcrXtalTds )
{
    //std::cout << "GcrXtalConvert::convert from ROOT to TDS BEGIN" << std::endl;   
    
    
    //CalXtalRecData* xTalRoot = gcrXtalRoot.getXtal();
    idents::CalXtalId tdsId = RootPersistence::convert(gcrXtalRoot.getXtalId()) ;    
    double pathLength = gcrXtalRoot.getPathLength();
    double closestFaceDist = gcrXtalRoot.getClosestFaceDist();
    int crossedFaces = gcrXtalRoot.getCrossedFaces();  
        
   // Event::CalXtalRecData* xTalTds = new Event::CalXtalRecData();
   // convert(*xTalRoot, *xTalTds);

    gcrXtalTds.setXtalId(tdsId);
    gcrXtalTds.setPathLength(pathLength);
    gcrXtalTds.setClosestFaceDist(closestFaceDist);
    gcrXtalTds.setCrossedFaces(crossedFaces);
    
    TVector3 entryPointRoot = gcrXtalRoot.getEntryPoint();
    TVector3 exitPointRoot = gcrXtalRoot.getExitPoint();

    Point    entryPointTds(entryPointRoot.x(), entryPointRoot.y(), entryPointRoot.z());
    Point    exitPointTds(exitPointRoot.x(), exitPointRoot.y(), exitPointRoot.z());

    gcrXtalTds.setEntryPoint(entryPointTds);
    gcrXtalTds.setExitPoint(exitPointTds);

    
    //std::cout << "GcrXtalConvert::convert from ROOT to TDS END" << std::endl;   
 
  
}

}
