// $Header$

/** @file MootData.cxx
    @author J. Bogart

    Non-inline implementation for moot-related objects in the tds.
*/

#include "CalibData/Moot/MootData.h"
#include "facilities/Util.h"

namespace CalibData {
  // copy constructors
  MootBase::MootBase(const MootBase& other)  {
    m_type = other.m_type;
    m_subtype = other.m_subtype;
  }
  
  MootParmCol::MootParmCol(MOOTSUBTYPE sub, unsigned key, unsigned sz)  :
    MootBase(MOOTTYPE_parm, sub), m_key(key) {
    
    if (sz) {
      m_v.reserve(sz);
      
    }
    
  }

  MootParmCol::MootParmCol(const MootParmCol& other) : MootBase(other), 
                                                       m_key(other.m_key)
  {
    m_v = other.m_v;
  }
  
  MootParmCol::~MootParmCol() {
    m_v.clear();
  }

  // MootFilterCfg
  MootFilterCfg::MootFilterCfg(const std::string& keyStr,
                               const std::string& name,
                               const std::string& pkg, 
                               const std::string& pkgVersion,
                               const std::string& fmxPath, 
                               const std::string& srcPath,
                               const std::string& fswIdStr,
                               const std::string& status,
                               const std::string& schemaIdStr,
                               const std::string& schemaVersionIdStr,
                               const std::string& instanceIdStr) :
    m_keyStr(keyStr), m_name(name), m_pkg(pkg), m_pkgVersion(pkgVersion),
    m_fmxPath(fmxPath), m_srcPath(srcPath), m_fswIdStr(fswIdStr), 
    m_status(status),
    m_schemaIdStr(schemaIdStr), m_schemaVersionIdStr(schemaVersionIdStr), 
    m_instanceIdStr(instanceIdStr) {
    
    using facilities::Util;
    
    m_key = Util::stringToUnsigned(m_keyStr);
    m_fswId = Util::stringToUnsigned(m_fswIdStr);
    m_schemaId = Util::stringToUnsigned(m_schemaIdStr);
    m_schemaVersionId = Util::stringToUnsigned(m_schemaVersionIdStr);
    m_instanceId = Util::stringToUnsigned(m_instanceIdStr);
  }
}     // end namespace CalibData
