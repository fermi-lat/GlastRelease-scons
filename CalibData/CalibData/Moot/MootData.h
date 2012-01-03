// $Header$

#ifndef CalibData_MootData_h
#define CalibData_MootData_h

#include <string>
#include "GaudiKernel/ContainedObject.h"
#include "GaudiKernel/ObjectVector.h"

/** @file MootData.h
    @author J. Bogart

    Definition of MOOT-related objects handled by MootSvc (however they
    don't go through Gaudi conversion service
*/

class MootSvc;


namespace CalibData {
  typedef enum {
    MOOTTYPE_notype,
    MOOTTYPE_parm,
    MOOTTYPE_anc,
    MOOTTYPE_vote
  } MOOTTYPE;
  
  typedef enum {
    MOOTSUBTYPE_nosubtype,
    MOOTSUBTYPE_latcParm,
    MOOTSUBTYPE_lciParm,
    MOOTSUBTYPE_cdmParm
  } MOOTSUBTYPE;
  
  class MootBase {
    friend class ::MootSvc;
  public:
    
    MootBase(MOOTTYPE t=MOOTTYPE_notype, MOOTSUBTYPE sub=MOOTSUBTYPE_nosubtype) :
      m_type(t), m_subtype(sub) {}
    MootBase(const MootBase& other);
    
    virtual ~MootBase() {};
    

    MOOTTYPE getType() const {return m_type;}
    MOOTSUBTYPE getSubtype() const {return m_subtype;}
    
  protected:
    MOOTTYPE m_type;
    MOOTSUBTYPE m_subtype;
    
  };  // end MootBase

  /**
     @class MootFilterCfg

     Pure data class to describe filter configuration cdm
   */
  class MootFilterCfg {
    friend class ::MootSvc;
  public:
    MootFilterCfg(const std::string& keyStr="", const std::string& name="",
                  const std::string& pkg="", const std::string& pkgVersion="",
                  const std::string& fmxPath="", const std::string& srcPath="",
                  const std::string& fswIdStr="", 
                  const std::string& status="",
                  const std::string& schemaIdStr="",
                  const std::string& schemaVersionIdStr="",
                  const std::string& instanceIdStr="");
    unsigned getKey() const {return m_key;}
    std::string getKeyStr() const {return m_keyStr;}
    std::string getName() const {return m_name;}
    std::string getPkg() const {return m_pkg;}
    std::string getPkgVersion() const {return m_pkgVersion;}
    std::string getFmxPath() const {return m_fmxPath;}
    std::string getSrcPath() const {return m_srcPath;}
    unsigned getFswId() const {return m_fswId;}
    std::string getFswIdStr() const {return m_fswIdStr;}
    std::string getStatus() const {return m_status;}
    unsigned getSchemaId() const {return m_schemaId;}
    std::string getSchemaIdStr() const {return m_schemaIdStr;}
    unsigned getSchemaVersionId() const {return m_schemaVersionId;}
    std::string getSchemaVersionIdStr() const {return m_schemaVersionIdStr;}
    unsigned getInstanceId() const {return m_instanceId;}
    std::string getInstanceIdStr() const {return m_instanceIdStr;}
  private:
    unsigned m_key;
    std::string m_keyStr;
    std::string m_name;
    std::string m_pkg;
    std::string m_pkgVersion;
    std::string m_fmxPath;
    std::string m_srcPath;
    unsigned m_fswId;
    std::string m_fswIdStr;
    std::string m_status;
    unsigned m_schemaId;
    std::string m_schemaIdStr;
    unsigned m_schemaVersionId;
    std::string m_schemaVersionIdStr;
    unsigned m_instanceId;
    std::string m_instanceIdStr;
  };
  
  /**
     @class MootParm
     
     Describes a single source file registered with moot as a "parameter"
  */

  class MootParm  {
    friend class ::MootSvc;
  public:

  
    MootParm(const std::string& key="", const std::string& pclass="", 
             const std::string& classFk="", const std::string& src="", 
             const std::string& srcFmt="", const std::string& status="",
             const std::string& precinct="") :
      m_key(key), m_class(pclass), m_classFk(classFk), m_src(src),
      m_srcFmt(srcFmt), m_status(status), m_precinct(precinct) {}
    
    std::string getKey() const {return m_key;}
    std::string getClass() const {return m_class;}
    std::string getClassFk() const {return m_classFk;}
    std::string getSrc() const {return m_src;}
    std::string getSrcFmt() const {return m_srcFmt;}
    std::string getStatus() const {return m_status;}
    std::string getPrecinct() const {return m_precinct;}
  
  private:
    std::string m_key;
    std::string m_class;
    std::string m_classFk;
    std::string m_src;
    std::string m_srcFmt;
    std::string m_status;
    std::string m_precinct;
  
  };      // end MootParm

  typedef std::vector<MootParm> MootParmVec;

  class MootParmCol : MootBase {

  friend class ::MootSvc;

  public:

    MootParmCol(MOOTSUBTYPE sub=MOOTSUBTYPE_nosubtype, unsigned key=0, 
                unsigned sz=0);

    MootParmCol(const MootParmCol& other);

    /// Client can use this to check whether information of interest might have
    /// changed recently. 
    unsigned fswKey() const {return m_key;}

    virtual ~MootParmCol();


    const MootParmVec& getMootParmVec() const {return m_v;}

  protected:
    MootParmVec m_v;
    unsigned   m_key;   // Exactly which key this is depends on MOOTSUBTYPE


  };
}
#endif
