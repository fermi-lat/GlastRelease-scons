// $Header$

#ifndef CalibData_MetadataInfo_H
#define CalibData_MetadataInfo_H

#include <iostream>
#include <string>
#include "CalibCnv/ICalibCnvSvc.h"

/** @class MetadataInfo

  Save as "real data" (in the TDDS) the stuff normally used to 
  look up the real data as a relatively simple test implementation
  of TDDS data + converter.

@author Joanne Bogart
 $Header$

*/

//extern const CLID& CLID_Calib_MetadataInfo;

namespace CalibData {
  class MetadataInfo : virtual public DataObject {

  public:

    MetadataInfo(const std::string& name, const std::string& value) :
      DataObject(), m_name(name), m_value(value) {}

    // Having these inline in include file could cause problems, in
    // fact the static member could already be a problem.  Will code
    // linked into different shareables have different copies of
    // 
    virtual const CLID& clID() const {return MetadataInfo::classID(); }
    static const CLID& classID();

    inline const std::string& getDataFmt() const {return m_dataFmt;}
    inline const std::string& getFmtVersion() const {return m_fmtVersion;}
    inline const std::string& getDataIdent() const {return m_dataIdent;}
    inline const ICalibCnvSvc::CalibKey getKey() const {return m_key;}

    virtual std::ostream& fillStream(std::ostream& s) const;

  private:
    ICalibCnvSvc::CalibKey   m_key;
    std::string m_dataFmt;
    std::string m_fmtVersion;
    std::string m_dataIdent;

  };
}





