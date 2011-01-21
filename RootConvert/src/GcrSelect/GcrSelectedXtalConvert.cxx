
#include <RootConvert/GcrSelect/GcrSelectedXtalConvert.h>
#include <RootConvert/Utilities/Toolkit.h>

namespace RootPersistence {

void convert( const Event::GcrSelectedXtal& gcrSelXtalTds, GcrSelectedXtal& gcrSelXtalRoot )
{

    CalXtalId rootId = RootPersistence::convert(gcrSelXtalTds.getXtalId()) ; 
    double pathLength =   gcrSelXtalTds.getPathLength ();
    double rawEnergy =   gcrSelXtalTds.getRawEnergy();
    double corrEnergy =   gcrSelXtalTds.getCorrEnergy ();
    double closestFaceDist = gcrSelXtalTds.getClosestFaceDist();
    int crossedFaces = gcrSelXtalTds.getCrossedFaces();  
    
    
    gcrSelXtalRoot.setXtalId(rootId);
    gcrSelXtalRoot.setPathLength(pathLength);
    gcrSelXtalRoot.setRawEnergy(rawEnergy);
    gcrSelXtalRoot.setCorrEnergy(corrEnergy);
    gcrSelXtalRoot.setSelectGrade(-1);
    
    gcrSelXtalRoot.setClosestFaceDist(closestFaceDist);
    gcrSelXtalRoot.setCrossedFaces(crossedFaces);

    Point entryPointTds = gcrSelXtalTds.getEntryPoint();
    Point exitPointTds = gcrSelXtalTds.getExitPoint();
    
    TVector3 entryPointRoot(entryPointTds.x(), entryPointTds.y(), entryPointTds.z());
    TVector3 exitPointRoot(exitPointTds.x(), exitPointTds.y(), exitPointTds.z());
    
    gcrSelXtalRoot.setEntryPoint(entryPointRoot);
    gcrSelXtalRoot.setExitPoint(exitPointRoot);

        /**std::cout << "gcrSelXtalTds->getXtalId()=" 
        << " @@TDS:@@@@@@"        
        << ", rawEnergy=" << gcrSelXtalTds.getRawEnergy()
        << ", corrEnergy=" << gcrSelXtalTds.getCorrEnergy()
        << ", crossedFaces=" << gcrSelXtalTds.getCrossedFaces()
        << ",entry.x()=" << gcrSelXtalTds.getEntryPoint().x()
        
        << " @@ROOT:@@@@@@"
        <<", rawEnergy=" << gcrSelXtalRoot.getRawEnergy()
        << ", corrEnergy=" << gcrSelXtalRoot.getCorrEnergy()
        << ", crossedFaces=" << gcrSelXtalRoot.getCrossedFaces() 
        << ",entry.x()=" << gcrSelXtalRoot.getEntryPoint().x()
        << endreq;*/

    
      //std::cout << "GcrXtalConvert::convert from TDS to ROOt END" << std::endl;   

}
 
void convert( const GcrSelectedXtal& gcrSelXtalRoot, Event::GcrSelectedXtal& gcrSelXtalTds )
{
    idents::CalXtalId tdsId = RootPersistence::convert(gcrSelXtalRoot.getXtalId()) ;    
    double pathLength = gcrSelXtalRoot.getPathLength();
    double rawEnergy =   gcrSelXtalRoot.getRawEnergy();
    double corrEnergy =   gcrSelXtalRoot.getCorrEnergy ();

    double closestFaceDist = gcrSelXtalRoot.getClosestFaceDist();
    int crossedFaces = gcrSelXtalRoot.getCrossedFaces();  
        
   // Event::CalXtalRecData* xTalTds = new Event::CalXtalRecData();
   // convert(*xTalRoot, *xTalTds);

    gcrSelXtalTds.setXtalId(tdsId);
    gcrSelXtalTds.setPathLength(pathLength);
    gcrSelXtalTds.setRawEnergy(rawEnergy);
    gcrSelXtalTds.setCorrEnergy(corrEnergy);
    gcrSelXtalTds.setSelectGrade(-1);
    
    gcrSelXtalTds.setClosestFaceDist(closestFaceDist);
    gcrSelXtalTds.setCrossedFaces(crossedFaces);
    
    TVector3 entryPointRoot = gcrSelXtalRoot.getEntryPoint();
    TVector3 exitPointRoot = gcrSelXtalRoot.getExitPoint();

    Point    entryPointTds(entryPointRoot.x(), entryPointRoot.y(), entryPointRoot.z());
    Point    exitPointTds(exitPointRoot.x(), exitPointRoot.y(), exitPointRoot.z());

    gcrSelXtalTds.setEntryPoint(entryPointTds);
    gcrSelXtalTds.setExitPoint(exitPointTds);

     
    //std::cout << "GcrXtalConvert::convert from ROOT to TDS END" << std::endl;   
 
  
}

}
