#include "CalibCLIDNode.h"
/**  @file CalibCLIDNode.cxx

     $Header$

     @author Joanne Bogart
*/

extern const CLID& CLID_Calib_CalibCLIDNode;
static const CLID& MetadataInfo::classID() {
    return CLID_Calib_CalibCLIDNode;
}

  std::ostream& MetadataInfo::fillStream(std::ostream& s) const {
    return s <<  "Child class ID = " << m_childClassID << std::endl;
  }
