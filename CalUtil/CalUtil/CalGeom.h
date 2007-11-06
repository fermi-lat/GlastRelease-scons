//    $Header$
#ifndef CalGeom_h
#define CalGeom_h

// LOCAL INCLUDES
#include "CalDefs.h"

// GLAST INCLUDES

// EXTLIB INCLUDES

// STD INCLUDES
#include <vector>

/** @file
    Geometry related utilities 
    @author Zach Fewtrell
*/

namespace CalUtil {
  /** generate list of active TEM modules in LAT */
  std::vector<CalUtil::TwrNum> findActiveTowers(IGlastDetSvc &detSvc);
};

#endif // CalGeom_h
