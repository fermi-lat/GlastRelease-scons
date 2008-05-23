// $Header$

#ifndef MOOT_FileDescrip_h
#define MOOT_FileDescrip_h

#include <string>

/** @file FileDescrip.h
    @author J. Bogart

*/
namespace MOOT {
  /**
        @class AcqSummaryInfo
   */
  class AcqSummaryInfo {
  public:
    AcqSummaryInfo(const std::string& key="",
                   const std::string& scid="",
                   const std::string& startedAt="",
                   const std::string& configKey="",
                   const std::string& acqType="",
                   const std::string& analysis="",
                   const std::string& status="",
                   const std::string& nevts="",
                   const std::string& evtutc0="",
                   const std::string& evtutc1="",
                   const std::string& hwKey="",
                   const std::string& swKey="",
                   const std::string& creator="",
                   const std::string& creationTime="",
                   const std::string& updateTime="") : 
      m_key(key), m_scid(scid), m_startedAt(startedAt), m_configKey(configKey),
      m_acqType(acqType), m_analysis(analysis), m_status(status), 
      m_nevts(nevts), m_evtutc0(evtutc0), m_evtutc1(evtutc1), m_hwKey(hwKey),
      m_swKey(swKey), m_creator(creator), m_creationTime(creationTime),
      m_updateTime(updateTime)
{}
    std::string getKey() const {return m_key;}
    std::string getScid() const {return m_scid;}
    std::string getStartedAt() const {return m_startedAt;}
    std::string getConfigKey() const {return m_configKey;}
    std::string getAcqType() const {return m_acqType;}
    std::string getAnalysis() const {return m_analysis;}
    std::string getStatus() const {return m_status;}
    std::string getNevts() const {return m_nevts;}
    std::string getEvtutc0() const {return m_evtutc0;}
    std::string getEvtutc1() const {return m_evtutc1;}
    std::string getHwKey() const {return m_hwKey;}
    std::string getSwKey() const {return m_swKey;}
    std::string getCreator() const {return m_creator;}
    std::string getCreationTime() const {return m_creationTime;}
    std::string getUpdateTime() const {return m_updateTime;}
  private:
    std::string m_key;
    std::string m_scid;
    std::string m_startedAt;
    std::string m_configKey;
    std::string m_acqType;
    std::string m_analysis;
    std::string m_status;
    std::string m_nevts;
    std::string m_evtutc0;
    std::string m_evtutc1;
    std::string m_hwKey;
    std::string m_swKey;
    std::string m_creator;
    std::string m_creationTime;
    std::string m_updateTime;
  };
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
               const std::string& verify="",
               const std::string& mode="", const std::string& cTime="",
               const std::string& voteKey="", const std::string& sbsKey="") : 
      m_key(key), m_name(name), m_alg(algorithm), m_step(step), 
      m_descrip(descrip), m_status(status), m_active(active),
      m_verify(verify),
      m_mode(mode), m_creationTime(cTime), m_voteKey(voteKey),
      m_sbsKey(sbsKey){};

    std::string getName() const {return m_name;}
    std::string getAlg() const {return m_alg;}
    std::string getStep() const {return m_step;}
    std::string getDescrip() const {return m_descrip;}
    std::string getStatus() const {return m_status;}
    std::string getActive() const {return m_active;}
    std::string getVerify() const {return m_verify;}
    std::string getKey() const {return m_key;}
    std::string getMode() const {return m_mode;}
    std::string getCreationTime() const {return m_creationTime;}
    std::string getVoteKey() const {return m_voteKey;}
    std::string getSbsKey() const {return m_sbsKey;}
    void clear() {
      m_key = m_name = m_alg = m_step = m_descrip = m_status = m_active 
        = m_verify = m_mode = m_creationTime = m_voteKey 
        = m_sbsKey = std::string("");
    }
  private:
    std::string m_key;
    std::string m_name;
    std::string m_alg;
    std::string m_step;
    std::string m_descrip;
    std::string m_status;
    std::string m_active;
    std::string m_verify;
    std::string m_mode;
    std::string m_creationTime;
    std::string m_voteKey;
    std::string m_sbsKey;
  };                                 // end ConfigInfo

  /** @class ConstitInfo
        Information concerning flight software constituents    
   */
  class ConstitInfo {
  public:
    ConstitInfo(const std::string& key = std::string(""),
                const std::string& name = std::string(""),
                const std::string& pkg = std::string(""),
                const std::string& pkgVersion = std::string(""),
                const std::string& dir = std::string(""),
                const std::string& fmxPath = std::string(""),
                const std::string& srcPath = std::string(""),
                const std::string& fswId = std::string(""),
                const std::string& status = std::string(""),
                const std::string& schemaId = std::string(""),
                const std::string& schemaVersionId = std::string(""),
                const std::string& instanceId = std::string("")) :
      m_key(key), m_name(name), m_pkg(pkg), m_pkgVersion(pkgVersion), 
      m_dir(dir), m_fmxPath(fmxPath), m_srcPath(srcPath), m_fswId(fswId), 
      m_status(status), m_schemaId(schemaId), 
      m_schemaVersionId(schemaVersionId), m_instanceId(instanceId)
    {};
    void clear() {
      m_key = m_name = m_pkgVersion = m_pkg = m_dir = m_fmxPath = m_srcPath
        = m_fswId = m_status = m_schemaId = m_instanceId
        = m_schemaVersionId = std::string("");
    }
    std::string getKey() const {return m_key;}
    std::string getName() const {return m_name;}
    std::string getPkg() const {return m_pkg;}
    std::string getPkgVersion() const {return m_pkgVersion;}
    std::string getDir() const {return m_dir;}
    std::string getFmxPath() const {return m_fmxPath;}
    std::string getSrcPath() const {return m_srcPath;}
    std::string getFswId() const {return m_fswId;}
    std::string getStatus() const {return m_status;}
    std::string getSchemaId() const {return m_schemaId;}
    std::string getSchemaVersionId() const {return m_schemaVersionId;}
    std::string getInstanceId() const {return m_instanceId;}
  private:
    std::string m_key;
    std::string m_name;
    std::string m_pkg;
    std::string m_pkgVersion;
    std::string m_dir;
    std::string m_fmxPath;
    std::string m_srcPath;
    std::string m_fswId;
    std::string m_status;
    std::string m_schemaId;
    std::string m_schemaVersionId;
    std::string m_instanceId;
  };

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
