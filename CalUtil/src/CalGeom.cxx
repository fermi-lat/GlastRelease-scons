/** @file
    implement CalGeom.h

    @author Zach Fewtrell
*/

// LOCAL INCLUDES
#include "CalUtil/CalDefs.h"

// GLAST INCLUDES
#include "idents/VolumeIdentifier.h"
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"

// EXTLIB INCLUDES
#include "GaudiKernel/StatusCode.h"

// STD INCLUDES
#include <vector>
#include <string>
#include <map>

namespace CalUtil {

  /// find which GLAST LAT towers are present by creating a VolumeIdentifer inside of each one and
  /// checking with GlastDetSvc if it is valid.
  std::vector<CalUtil::TwrNum>  findActiveTowers(IGlastDetSvc &detSvc) {
    std::vector<CalUtil::TwrNum> twrList;
  
    /*-- RETRIEVE CONSTANTS */
    double tmp;
    // map containing pointers to integer constants to be read
    // with their symbolic names from xml file used as a key 
    typedef std::map<int*,string> PARAMAP;
    PARAMAP param;

    int eTowerCAL, eLATTowers, nCsISeg, eXtal, eMeasureX;

    // INT CONSTANTS
    //     filling the map with information on constants to be read 
    param[&eTowerCAL]    = std::string("eTowerCAL");
    param[&eLATTowers]   = std::string("eLATTowers");
    param[&nCsISeg]      = std::string("nCsISeg");
    param[&eXtal]        = std::string("eXtal");
    param[&eMeasureX]    = std::string("eMeasureX");

    // loop over all constants information contained in the map
    for(PARAMAP::iterator iter=param.begin(); iter!=param.end();iter++){
      //  retrieve constant
      if(!detSvc.getNumericConstByName((*iter).second, &tmp))
        return twrList;
      else *((*iter).first)= int(tmp); // store retrieved value 
    }


    for (CalUtil::TwrNum testTwr; testTwr.isValid(); testTwr++) {
      // create geometry ID for 1st xtal in each tower
      idents::VolumeIdentifier volId;

      // volId info snagged from 
      // http://www.slac.stanford.edu/exp/glast/ground/software/geometry/docs/identifiers/geoId-RitzId.shtml
      volId.append(eLATTowers);
      volId.append(testTwr.getRow());
      volId.append(testTwr.getCol());
      volId.append(eTowerCAL);
      volId.append(0); // layer
      volId.append(eMeasureX);
      volId.append(0); // column
      volId.append(eXtal);
      volId.append(0); // segment Id

      // test to see if the volume ID is valid.
      std::string tmpStr; std::vector<double> tmpVec;
      StatusCode sc = detSvc.getShapeByID(volId, &tmpStr, &tmpVec);
      if (!sc.isFailure())
        twrList.push_back(testTwr);
    }
    
    return twrList;
  }
}
