#include "CalibData/MetadataInfo.h"
/**  @file MetadataInfo.cxx

     $Header$

     @author Joanne Bogart
*/

extern const CLID& CLID_Calib_MetadataInfo;
namespace CalibData {
  static const CLID& MetadataInfo::classID() {
    return CLID_MetadataInfo;
  }

  std::ostream& MetadataInfo::fillStream(std::ostream& s) const {
    return s <<  "Key = " << key << " data format: " << m_dataFmt 
             << " version: " << m_fmtVersion << " data ident: "
             << m_dataIdent << std::endl;
  }
