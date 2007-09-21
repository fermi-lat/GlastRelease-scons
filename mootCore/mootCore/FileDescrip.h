// $Header$

#ifndef MOOT_FileDescrip_h
#define MOOT_FileDescrip_h

#include <string>

/** @file FileDescrip.h
    @author J. Bogart

*/
namespace MOOT {
  /**
        @class AncAliasInfo
  */
  class AncAliasInfo {
  public:
    AncAliasInfo(const std::string& key="", const std::string& name="",
                 const std::string& aclass="", const std::string& aclassFk="", 
                 const std::string& ancFk="",
                 const std::string& lastModifiedTs="", 
                 const std::string& lastModifier="") :
      m_key(key), m_name(name), m_class(aclass), m_classFk(aclassFk), 
      m_ancFk(ancFk), m_rowTs(lastModifiedTs), m_modifier(lastModifier) {}

    std::string getKey() const {return m_key;}
    std::string getAncKey() const {return m_ancFk;}
    std::string getAliasName() const {return m_name;}
    std::string getClass() const {return m_class;}
    std::string getClassFk() const {return m_classFk;}
    std::string getRowModTs() const {return m_rowTs;}
    std::string getModifier() const {return m_modifier;}

  private:
    std::string m_key;
    std::string m_name;
    std::string m_class;
    std::string m_classFk;
    std::string m_ancFk;
    std::string m_rowTs;
    std::string m_modifier;
  };  // end AncAliasInfo

  /**
        @class AncInfo
  */
  class AncInfo {
  public:
    AncInfo(const std::string& key="", const std::string& aclass="", 
            const std::string& classFk="", const std::string& instr="",
            const std::string& src="", const std::string& srcFmt="",
            const std::string& fileCreationTs="", 
            const std::string& rowCreationTs="", 
            const std::string& rowCreator="", const std::string& quality="",
            const std::string& descrip="", const std::string& status="") :
      m_key(key), m_class(aclass), m_classFk(classFk), m_instr(instr),
      m_src(src), m_srcFmt(srcFmt), m_fileTs(fileCreationTs), 
      m_rowTs(rowCreationTs), m_creator(rowCreator), m_quality(quality),
      m_descrip(descrip), m_status(status) {}

    std::string getKey() const {return m_key;}
    std::string getClass() const {return m_class;}
    std::string getClassFk() const {return m_classFk;}
    std::string getInstr() const {return m_instr;}
    std::string getSrc() const {return m_src;}
    std::string getSrcFmt() const {return m_srcFmt;}
    std::string getFileCreation() const {return m_fileTs;}
    std::string getRowCreation() const {return m_rowTs;}
    std::string getCreator() const {return m_creator;}
    std::string getQuality() const {return m_quality;}
    std::string getDescrip() const {return m_descrip;}
    std::string getStatus() const {return m_status;}

  private:
    std::string m_key;
    std::string m_class;
    std::string m_classFk;
    std::string m_instr;
    std::string m_src;
    std::string m_srcFmt;
    std::string m_fileTs;
    std::string m_rowTs;
    std::string m_creator;
    std::string m_quality;
    std::string m_descrip;
    std::string m_status;
  };

  /**
        @class ConfigInfo
  */
  class ConfigInfo {
  public:
    ConfigInfo(const std::string& key="", const std::string& name="", 
               const std::string& algorithm="", const std::string& step="", 
               const std::string& descrip="",const std::string& status="",
               const std::string& active="",
               const std::string& mode="", const std::string& cTime="") : 
      m_key(key), m_name(name), m_alg(algorithm), m_step(step), 
      m_descrip(descrip), m_status(status), m_active(active),
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
  };                                 // end ConfigInfo

  /**
    Little utility class used by MootBuild and MootQuery.
    Define a type consisting of two strings, one for file path
    and one a string for type of file (e.g., "latc_CCC")
   */
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
    /// Currently (April, 2007) values in active use are of the form 
    ///  lci_<subsystem>  (i.e., lci_ACD, lci_CAL or lci_TKR)
    ///  ltc_C, ltc_F, ltc_T   (for thermal control)
    ///  latc_ignore
    ///  latc_<component>_backup 
    /// or latc_<component>_<precinct>  
    /// where <component> is one of the latc component names
    ///    DFT,GEM,AEM,ARC,AFE,TEM,TIC,CCC,CRC,CFE,TCC,TRC,SPT,TFE,TDC
    //  and <precinct> is a value in Precincts.name
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

  /**
       @class ParmInfo
   */
  class ParmInfo {
  public:
    ParmInfo(const std::string& key="", const std::string& pclass="", 
             const std::string& classFk="", const std::string& voteFk="",
             const std::string& instr="", const std::string& src="", 
             const std::string& srcFmt="", const std::string& rowTs="", 
             const std::string& rowCreator="", const std::string& quality="",
             const std::string& descrip="", const std::string& status="") :
      m_key(key), m_class(pclass), m_classFk(classFk), m_voteFk(voteFk),
      m_instr(instr), m_src(src), m_srcFmt(srcFmt), m_rowTs(rowTs), 
      m_creator(rowCreator), m_quality(quality), m_descrip(descrip), 
      m_status(status) {}

    std::string getKey() const {return m_key;}
    std::string getClass() const {return m_class;}
    std::string getClassFk() const {return m_classFk;}
    std::string getVoteFk() const {return m_voteFk;}
    std::string getInstr() const {return m_instr;}
    std::string getSrc() const {return m_src;}
    std::string getSrcFmt() const {return m_srcFmt;}
    std::string getRowCreation() const {return m_rowTs;}
    std::string getCreator() const {return m_creator;}
    std::string getQuality() const {return m_quality;}
    std::string getDescrip() const {return m_descrip;}
    std::string getStatus() const {return m_status;}

  private:
    std::string m_key;
    std::string m_class;
    std::string m_classFk;
    std::string m_voteFk;
    std::string m_instr;
    std::string m_src;
    std::string m_srcFmt;
    std::string m_rowTs;
    std::string m_creator;
    std::string m_quality;
    std::string m_descrip;
    std::string m_status;
  };             // end ParmInfo

  /**
    @class ParmOffline
     
    Information about parameter files of interest to Offline
  */
  class ParmOffline {
  public:
    ParmOffline(const std::string& key="", const std::string& pclass="", 
                const std::string& classFk="", const std::string& src="", 
                const std::string& srcFmt="", const std::string& status="") :
      m_key(key), m_class(pclass), m_classFk(classFk),
      m_src(src), m_srcFmt(srcFmt), m_status(status) {}

    std::string getKey() const {return m_key;}
    std::string getClass() const {return m_class;}
    std::string getClassFk() const {return m_classFk;}
    std::string getSrc() const {return m_src;}
    std::string getSrcFmt() const {return m_srcFmt;}
    std::string getStatus() const {return m_status;}

  private:
    std::string m_key;
    std::string m_class;
    std::string m_classFk;
    std::string m_src;
    std::string m_srcFmt;
    std::string m_status;

  };      // end ParmOffline
  /**
        @class voteAliasInfo
  */
  class VoteAliasInfo {
  public:
    VoteAliasInfo(const std::string& key="", const std::string& name="",
                  const std::string& precinct="", 
                  const std::string& precinctFk="", 
                  const std::string& voteFk="",
                  const std::string& lastModifiedTs="", 
                  const std::string& lastModifier="") :
      m_key(key), m_name(name), 
      m_precinct(precinct), m_precinctFk(precinctFk), m_voteFk(voteFk), 
      m_rowTs(lastModifiedTs), m_modifier(lastModifier) {}

    std::string getKey() const {return m_key;}
    std::string getAliasName() const {return m_name;}
    std::string getPrecinct() const {return m_precinct;}
    std::string getPrecinctFk() const {return m_precinctFk;}
    std::string getVoteKey() const {return m_voteFk;}
    std::string getRowModTs() const {return m_rowTs;}
    std::string getModifier() const {return m_modifier;}

  private:
    std::string m_key;
    std::string m_name;
    std::string m_precinct;
    std::string m_precinctFk;
    std::string m_voteFk;
    std::string m_rowTs;
    std::string m_modifier;
  };  // end VoteAliasInfo


  /**
        @class VoteInfo
  */
  class VoteInfo {
  public:
    VoteInfo(const std::string& key="", const std::string& precinct="",
             const std::string& precinctFk="", const std::string& instr="", 
             const std::string& src="", const std::string& ts="", 
             const std::string& creator="",
             const std::string& descrip="", const std::string& status="") :
      m_key(key), m_precinct(precinct), m_precinctFk(precinctFk), 
      m_instr(instr),
      m_src(src), m_ts(ts), m_creator(creator),
      m_descrip(descrip), m_status(status)  {};
    std::string getKey() const {return m_key;}
    std::string getPrecinctFk() const {return m_precinctFk;}
    std::string getPrecinct() const {return m_precinct;}
    std::string getSrc() const {return m_src;}
    std::string getStatus() const {return m_status;}
    std::string getInstr() const {return m_instr;}
    std::string getDescrip() const {return m_descrip;}
    std::string getTs() const {return m_ts;}
    std::string getCreator() const {return m_creator;}
      private:
    std::string m_key;
    std::string m_precinct;     // name
    std::string m_precinctFk;   // foreign key as string
    std::string m_instr;
    std::string m_src;
    std::string m_ts;
    std::string m_creator;
    std::string m_descrip;
    std::string m_status;
  };    // end VoteInfo
             

}

#endif
