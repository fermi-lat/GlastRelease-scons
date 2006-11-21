// $Header$

#ifndef MOOT_FileDescrip_h
#define MOOT_FileDescrip_h

#include <string>

/** @file FileDescrip.h
    @author J. Bogart

    Little utility class used by MootBuild and MootQuery.
    Define a type consisting of two strings, one for file path
    and one a string for type of file (e.g., "latc_CCC")

*/
namespace MOOT {
  class FileDescrip {
  public:
    FileDescrip(const std::string& path="", 
                const std::string& type="") : m_path(path), m_type(type) {}
    std::string getPath() const { return m_path;}
    std::string getType() const { return m_type;}
    void setPath(std::string path) {m_path = path;}
    void setType(std::string type) {m_type = type;}
    /// Where to find the parameter data; typically a file path
    std::string m_path;
    /// Parameter type; should match something in Parameter_class table
    /// Currently (April, 2006) allowable values are of the form 
    ///  lci_<subsystem>  (i.e., lci_ACD, lci_CAL or lci_TKR)
    /// or latc_<component>  
    /// where <component> is one of the latc component names
    ///    DFT,GEM,AEM,ARC,AFE,TEM,TIC,CCC,CRC,CFE,TCC,TRC,SPT,TFE,TDC
    std::string m_type;  
  };
  class FswInfo {
  public:
    FswInfo(const std::string& fmxPath="", const std::string& src="",
            const std::string& description="", 
            const std::string& fmxDb="", unsigned fmxKey=0,
            unsigned mootClassKey=0) : 
      m_fmxPath(fmxPath), m_src(src), m_desc(description), m_fmxDb(fmxDb),
      m_fmxKey(fmxKey), m_classKey(mootClassKey) {};
    std::string getFmxPath() const {return m_fmxPath;}
    std::string getSrc() const {return m_src;}
    std::string getDescription() const {return m_desc;}
    std::string getFmxDb() const {return m_fmxDb;}
    unsigned getFmxKey() const {return m_fmxKey;}
    unsigned getMootClassKey() const {return m_classKey;}
  private:
    std::string m_fmxPath;
    std::string m_src;
    std::string m_desc;
    std::string m_fmxDb;
    unsigned m_fmxKey;
    unsigned m_classKey;
  };


  class ConfigInfo {
  public:
    ConfigInfo(const std::string& key="", const std::string& name="", 
               const std::string& algorithm="", const std::string& step="", 
               const std::string& description="",const std::string& status="",
               const std::string& active="",
               const std::string& mode="", const std::string& cTime="") : 
      m_key(key), m_name(name), m_alg(algorithm), m_step(step), 
      m_descrip(description), m_status(status), m_active(active),
      m_mode(mode), m_creationTime(cTime) {};

    std::string getName() const {return m_name;}
    std::string getAlg() const {return m_alg;}
    std::string getStep() const {return m_step;}
    std::string getDescrip() const {return m_descrip;}
    std::string getStatus() const {return m_status;}
    std::string getActive() const {return m_active;}
    std::string getKey() const {return m_key;}
    std::string getMode() const {return m_mode;}
    std::string getCreationTime() const {return m_creationTime;}
  private:
    std::string m_key;
    std::string m_name;
    std::string m_alg;
    std::string m_step;
    std::string m_descrip;
    std::string m_status;
    std::string m_active;
    std::string m_mode;
    std::string m_creationTime;
  };

  /*
     @class InputDescrip
 
     Similar to FileDescrip, but here include extra enum field to
     distinguish sort of input: parameter file, ancillary file or 
     offline calibraton.
  */
  /*
  class InputDescrip {         don't know if we need this yet
  public:
    InputDescrip(REGISTERTYPE       what,
                 const std::string& path,
                 const std::string& type)
      :  m_registerType(what), m_path(path), m_type(type) {}

    /// What kind of thing is being registered
    REGISTERTYPE m_registerType;
    /// Where to find the data; typically a file path
    std::string m_path;
    /// Parameter type; should match something in Parameter_class table
    std::string m_type;  
  };
  */

}

#endif
