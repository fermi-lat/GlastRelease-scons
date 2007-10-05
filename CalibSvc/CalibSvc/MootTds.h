// $Header$

#ifndef MootTds_h
#define MootTds_h

#include <string>
#include "GaudiKernel/ContainedObject.h"
#include "GaudiKernel/ObjectVector.h"

/** @file MootTds.h
    @author J. Bogart

    Definition of MOOT-related objects stored in the "calibration" TDS
*/

static const CLID CLID_MootBase = 10000;
static const CLID CLID_Moot = CLID_MootBase;
static const CLID CLID_MootParm = CLID_Moot + 1;


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

class MootBase : public DataObject {
  friend class MootBaseCnv;  

public:
  
  MootBase(MOOTTYPE t=MOOTTYPE_notype, MOOTSUBTYPE sub=MOOTSUBTYPE_nosubtype) :
    m_type(t), m_subtype(sub) {}
  MootBase(const MootBase& other);
  
  virtual ~MootBase() {};

  // Re-implemented from DataObject
  /// Class ID of this instance
  inline virtual const CLID& clID() const { return classID(); } 
  
  /// Class ID of this class
  inline static  const CLID& classID() { return CLID_MootBase; };

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

class MootParm : virtual public ContainedObject {
public:
  virtual const CLID& clID() const   { return MootParm::classID(); }
  static const CLID& classID()       { return CLID_MootParm; }
  
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


typedef ObjectVector<MootParm> MootParmOV;

class MootParmCol : public MootBase {
  friend class MootParmCnv;
public:

  MootParmCol(MOOTSUBTYPE sub=MOOTSUBTYPE_nosubtype, unsigned key=0, 
              unsigned sz=0);

  MootParmCol(const MootParmCol& other);

  /// Client can use this to check whether information of interest might have
  /// changed recently. 
  unsigned fswKey() const {return m_key;}

  virtual ~MootParmCol();

  // Re-implemented from DataObject
  inline virtual const CLID& clID() const { return classID(); } 
  inline static  const CLID& classID() { return CLID_MootParm; };

  const MootParmOV& getMootParmOV() const {return m_ov;}

protected:
  MootParmOV m_ov;
  unsigned   m_key;   // Exactly which key this is depends on MOOTSUBTYPE


};

#endif
