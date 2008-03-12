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
    friend class MootSvc;
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
     @class MootParm
     
     Describes a single source file registered with moot as a "parameter"
  */

  class MootParm  {
    friend class MootSvc;
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

  class MootParmCol : public MootBase {
    friend class MootSvc;
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
