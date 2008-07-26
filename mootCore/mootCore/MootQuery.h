// $Header$
// Handles registering a single file into Parameters table.  
// Might also handle other parts of building a config; TBD.
#ifndef MOOT_MootQuery_h
#define MOOT_MootQuery_h

#include <string>
#include <vector>
#include <set>

#include "rdbModel/Tables/Column.h"
#include "mootCore/FileDescrip.h"

#include <xercesc/util/XercesDefs.hpp>

#ifndef SWIG
XERCES_CPP_NAMESPACE_BEGIN
#endif
class  DOMElement;

#ifndef SWIG
XERCES_CPP_NAMESPACE_END
#endif

/** @file MootQuery.h
    @author J. Bogart
  
    Declare the class MootQuery which provides facilities to 
    register components of a configuration in the configuration tracking
    database MOOD.
*/
namespace rdbModel {
  class Rdb;
  class Connection;
}

namespace facilities {
  class Timestamp;
}

typedef std::vector<std::string> VectorOfString;
namespace MOOT {
  class MoodConnection;
  class FilterAssoc; // Def is in filterConfig.cxx, only place it's used

  enum LpaMode {
    LPA_MODE_NORMAL = 0,
      LPA_MODE_TOO = 1,
    LPA_MODE_GRB0 = 2,
    LPA_MODE_GRB1 = 3,
    LPA_MODE_GRB2 = 4,
    LPA_MODE_SOLAR = 5,
    LPA_MODE_CALIBRATION = 6,
    LPA_MODE_DIAGNOSTIC = 7,
    LPA_MODE_count,
    LPA_MODE_ALL = 0x8000 
  };

  // assigned by MOC.  This is default for getAcqSummaryInfo.
  // For MBR use 80
  const int FLIGHT_SCID = 77; 

  class MootQuery   {
  public:
    /// Use connection if one is supplied; else get a new read connection
    MootQuery(MoodConnection* mood=0);

    ~MootQuery();

    /**
       Translate string to key.  Return -1 if there is no translation
     */
    int classKey(const std::string& table, const std::string& name);

    /**
       Translate  class key to corresponding string.
       Return empty string if there is no translation
    */
    std::string classStr(const std::string& table, unsigned key);


    int fswClassKey(const std::string& name) {
      return classKey("FSW_class", name);
    }
    std::string fswClassString(unsigned key) {
      return classStr("FSW_class", key);
    }

    /** Store everything of interest about this row in Acq_summary
        table in AcqSummaryInfo structure.
        @return 0 if no row exists with this key;  otherwise ptr
         to filled structure.  Caller is responsible for deleting.
     */
    AcqSummaryInfo* getAcqSummaryInfoByKey(unsigned key);

    /**
       (scid, startedAt) is another way to uniquely identify an acquisition
     */
    AcqSummaryInfo* getAcqSummaryInfo(unsigned startedAt, 
                                      unsigned scid=FLIGHT_SCID);

    /**
       Find acqu summary specified by an event time and scid
     */
    AcqSummaryInfo* getAcqSummaryInfo(const facilities::Timestamp& eventTime,
                                      unsigned scid=FLIGHT_SCID);


    /** Store everything of interest about this row in AncillaryAlias
        table in AncAliasInfo structure.
        @return 0 if no row exists with this key;  otherwise ptr
         to filled structure.  Caller is responsible for deleting.
     */
    AncAliasInfo* getAncAliasInfo(unsigned key);


    /**
       Return ConstitInfo objects for all active filters in specified
       mode.  Return union of active filters for all modes for special
       value LPA_MODE_ALL (the default)
       By default the input vector will be cleared first.
     */
    unsigned getActiveFilters(unsigned configKey, 
                              std::vector<ConstitInfo>& info,
                              LpaMode mode = LPA_MODE_ALL,
                              bool clear=true);

    /**
       Return pointer to ConstitInfo object, given configKey,
       lpa mode and handler id if there is a filter config satisfying
       these constraints; else return null pointer.
       @a handlerName  is also returned (empty string if none)
     */
    const ConstitInfo* getActiveFilter(unsigned configKey, 
                                       LpaMode mode,
                                       unsigned handlerId,
                                       std::string& handlerName);

    unsigned getAncClasses(std::vector<std::string>& names);

    /** Store everything of interest about this row in Ancillary
        table in AncInfo structure.
        @return 0 if no row exists with this key; otherwise ptr
         to filled structure.  Caller is responsible for deleting.
     */
    AncInfo* getAncInfo(unsigned key);

    /** Return list of keys of anc. files associated with one or more
        of the input parameter file keys.  Return count (may include
        duplicates).
     */
    unsigned getAncsFromParms(const std::vector<unsigned>& parms,
                              std::vector<unsigned>& ancs);

    /**
       Return keys of ancillary files which were input to parameter entries 
       supplied when config was created.
     */
    bool getConfigAncsRequest(unsigned configKey,
                              std::vector<unsigned>& ancKeys);
   /**
      Return keys of ancillary files which were input to parameter entries
      which were in turn input for binaries belonging to the config.
   */
    bool getConfigAncsUsed(unsigned configKey,
                           std::vector<unsigned>& ancKeys);
    /**
       Given a config key, return associated list of fmx relative
       file paths.  Such a file path is what fmx upload expects
       as input.  Require Config to have good status.
     */
    bool getConfigFmxPaths(const std::string& configKey, 
                           std::vector<std::string>& paths);

    bool getConfigFmxPaths(unsigned configKey, 
                           std::vector<std::string>& paths);
    /**
       look up by alg name and step; use most recent config with good status
     */
    bool getConfigFmxPathsByAlg(const std::string& algName,
                                const std::string& algStep,
                                std::vector<std::string>& paths);

    bool getConfigFmxPathsByAlg(const std::string& algName,
                                unsigned algStep,
                                std::vector<std::string>& paths);
    /**
       look up by name; use most recent config with good status
     */
    bool getConfigFmxPathsByName(const std::string& configName, 
                                 std::vector<std::string>& paths);
    bool getConfigFswInfo(const std::string& configKey,
                          std::vector<FswInfo>& info, bool clear=true);


    bool getConfigFswInfo(unsigned configKey,
                          std::vector<FswInfo>& info, bool clear=true);

    bool getConfigFswInfoByAlg(const std::string& algName,
                               const std::string& algStep,
                                std::vector<FswInfo>& info, bool clear=true);


    bool getConfigFswInfoByAlg(const std::string& algName,
                               unsigned algStep,
                               std::vector<FswInfo>& info, bool clear=true);

    bool getConfigFswInfoByName(const std::string& configName,
                                std::vector<FswInfo>& info, bool clear=true);
    /**
       Store a list of structures in @info argument concerning configs 
       which satisfy cuts on status, instrument and mode.  They're sorted
       in ascending order by key.  Return value is count of Configs
       satisfying the cuts.
     */
    unsigned getConfigInfo(std::vector<ConfigInfo>& info,
                           const std::string& status="CREATED",
                           const std::string& instr="LAT",
                           const std::string& mode="");

    /**
       Return pointer to ConfigInfo for the single Config specified by key.
       Return 0 if not found or other failure.
     */
    ConfigInfo* getConfigInfo(unsigned key);
    /**
       Given a config key (unsigned), return associated list of fsw input
       keys (that is, primary keys from MOOT's FSW_inputs table)
       By default, only return list if the Config has status
       "CREATED".
     */
    bool getConfigInputs(unsigned configKey, 
                         std::vector<std::string>& inputKeys,
                         bool goodStatus = true);

    /**
       Given a config key (string), return associated list of fsw input
       keys (that is, primary keys from MOOT's FSW_inputs table)
       By default, only return list if the Config has status
       "CREATED".
     */
    bool getConfigInputs(const std::string& configKeyStr, 
                         std::vector<std::string>& inputKeys,
                         bool goodStatus = true);

    /**
       return list of all keys of configs with algorithm name algname,
       alg_step = step.  Empty string for @a status or @a instr 
       (or @a name argument for getConfigKeysByName) means any value 
       is acceptable. Similarly, @a step=0 means Config can have any 
       value for step.
     */
    unsigned getConfigKeysByAlg(std::vector<unsigned>& keys, 
                                const std::string& algname,
                                unsigned step=1,
                                const std::string& status="CREATED",
                                const std::string& instr="LAT");

    unsigned getConfigKeysByName(std::vector<unsigned>& keys, 
                                 const std::string& name, 
                                 const std::string& status="CREATED",
                                 const std::string& instr="LAT");
    /**
       Find keys of parameter entries supplied when config was
       created.
     */
    bool getConfigParmsRequest(unsigned configKey, 
                               std::vector<unsigned>& parameterKeys);

    /**
       Find keys of parameter entries used to build FSW inputs
       for this config
     */
    bool getConfigParmsUsed(unsigned configKey, 
                            std::vector<unsigned>& parameterKeys);

    /**
       Find list of keys for all Configs using vote with alias
       @a alias for precinct @a precinct.
     */
    bool getConfigsForAlias(const std::string& precinct, 
                            const std::string& alias,
                            std::set<unsigned>& keys);

    /**
       Return list of keys for all Configs with Parameter files dependent
       on the ancillary register with key @a ancKey
    */
    bool getConfigsForAncillary(unsigned ancKey, 
                                std::set<unsigned>& keys);


    /**
       Return list of keys for all Configs whose Parameter files were
       built using specified vote
    */
    bool getConfigsForVote(unsigned voteKey, std::set<unsigned>& keys);


    ConstitInfo* getConstituentInfo(unsigned constitKey);

    ConstitInfo* getConstituentByFswId(unsigned fswId);

    /**
      Return info for all constituents associated with the 
      supplied @a sbsFmxExportKey. If @a schemaId is non-negative,
      only return info for constituents with that value for schema_id.
      @return  count of constituents satisfying constraints
    */
    unsigned getConstituentInfoBySbs(unsigned sbsFmxExportKey, 
                                     std::vector<ConstitInfo>& info,
                                     int schemaId=-1);


    unsigned getLastConfigKeyByAlg(const std::string& alg, unsigned step=1,
                                   const std::string& status="CREATED",
                                   const std::string& instr="LAT");


    unsigned getLastConfigKeyByName(const std::string& name, 
                                    const std::string& status="CREATED",
                                    const std::string& instr="LAT");


    /**
       Given fmx log. key for latc master file, return list of
       pairs. Each element in pair consists of file path in
       MOOT archive and string specifying latc component
       ..what about tower #?  Might need that if LATC inputs
        are split along tower boundaries
    */
    bool getLatcSrc(const std::string& latcMasterKey,
                    std::vector<FileDescrip>& sources);

    bool getLatcSrc(unsigned latcMasterKey,
                    std::vector<FileDescrip>& sources);

    /** 
        Given SBS logical key, return constituent info for all
        lpa_db constituents
     */
    bool getLpaConstituents(unsigned sbsKey, std::vector<ConstitInfo>& lpas);


    /**
       Get LATC master (i.e., fmx logical key for it) belonging to
       @a configKey
     */
    unsigned getMasterKey(unsigned configKey);
    /** 
        Given SBS logical key, return constituent info for all
        constituents in the specified package
     */
    bool getPackageConstituents(unsigned sbsKey, const std::string& pkg,
                                std::vector<ConstitInfo>& constits);

    /** Store everything of interest about this row in Parameters
        table in Info structure.
        @return 0 if no row exists with this key; 1 otherwise
     */
    ParmInfo* getParmInfo(unsigned key);


    /// Fill supplied argument with parameter class names. Return count
    unsigned getParmClasses(std::vector<std::string>& names, 
                            const std::string& precinct="*");

    /// Given fmx logical id for a latc master, look up associated
    /// parameter files (latc source) and return a little structure
    /// of information about each.
    bool getParmsFromMaster(unsigned fmxMasterKey, 
                            std::vector<ParmOffline>& parms);


    /// Old name for getParmClassees
    unsigned getParameterClasses(std::vector<std::string>& names);

    unsigned getPrecincts(std::vector<std::string>& names);

    VoteAliasInfo* getVoteAliasInfo(unsigned key);

    VoteInfo* getVoteInfo(unsigned key);

    /** If vote file specified by @a voteKey is up to date, put keys 
       of parameter files of all relevant classes associated with the
       vote file in the @a parmKeys vector and return true.  
       Else 
            if (!force) return false and an empty list. 
            if (force) {
                  still ultimately return false, but don't clear list
                  [for container, if any contained vote has no
                   parameters at all, clear list]
            }
              
     */
    bool getVoteParameters(unsigned voteKey, 
                           std::vector<unsigned>& parmKeys, bool force=false);

    /**
       Get all info for registered votes belonging to specified
       precinct
     */
    bool getVotesForPrecinct(const std::string& precinct, 
                             std::vector<VoteInfo>& voteInfo);

    /**
       Return ancillary alias keys satisfying specified conditions on fields
       in Ancillary_alias table.  "*" means "don't cut on this condition".
       Keys are returned in descending order.
     */
    unsigned listAncAliasKeys(std::vector<unsigned>& keys,
                              const std::string& aClass="*",
                              const std::string& aliasName="*");

    /**
       Return ancillary keys satisfying specified conditions on fields
       in Ancillary table.  "*" means "don't cut on this condition"
     */
    unsigned listAncKeys(std::vector<unsigned>& keys,
                         const std::string& aClass="*",
                         const std::string status="CREATED",
                         const std::string quality="PROD",
                         const std::string& instr="LAT");

    /**
       Return list of keys for all good Configs (by default) or all
       configs with specified status.  By convention, if status
       argument is supplied and is equal to empty string, all
       config keys will be returned.
     */
    unsigned listConfigKeys(std::vector<unsigned>& keys, 
                            const std::string& status="CREATED",
                            const std::string& instr="LAT");
                  
    /**
       Return vote alias keys satisfying specified conditions on fields
       in Vote_aliases table.  "*" means "don't cut on this condition"
     */
    unsigned listVoteAliasKeys(std::vector<unsigned>& keys,
                               const std::string& precinct="*",
                               const std::string& aliasName="*");
          
    /**
       Return vote keys satisfying specified conditions on fields
       in Votes table.  "*" means "don't cut on this condition"
     */
    unsigned listVoteKeys(std::vector<unsigned>& keys,
                          const std::string& precinct="*",
                          const std::string& status="CREATED",
                          const std::string& instr="LAT");


    /// Returns key of file registered in Ancillary if there is one; else 0.
    /// @a ancClass is ancillary class name.
    unsigned resolveAncAlias(const std::string& alias, 
                             const std::string& ancClass);

    /** 
        Return number of ancillary file keys found.  Throw exception
        for error (e.g., no such voteKey; no ancillary found for an alias) 
     */
    bool resolveAncAliases(std::vector<std::string>& ancKeys,
                           unsigned voteKey);

    /// Returns key of file registered in Votes if there is one; else 0
    unsigned resolveVoteAlias(const std::string& alias,
                              const std::string& precinct);

    /**
       Return false if
          - voteKey doesn't exist in db
          - voteKey does exist, but @arg goodStatus was true and 
            vote did not have status "CREATED"
    */
    bool voteExists(unsigned voteKey, bool goodStatus=true);

    /// The "vote" referred to may concern a single precinct or may
    /// be container file, referring to several individual precinct files.
    /// If optional arg. is supplied and vote is up to date return
    /// parm keys associated with the vote.
    /// If isCtn is non-zero, inform caller whether or not vote was container

    bool voteIsUpToDate(unsigned voteKey, std::vector<unsigned>* parmKeys=0,
                        bool* isCtn=0, bool force=false);

  private:
    bool cacheFilterConfig(unsigned configKey);
    bool parseFilterParm(const std::string& path);
    /// Parse mode element in a filter assoc.  parameter file and
    /// update cache accordingly
    bool parseLpaModeElt(
#ifndef SWIG
                         XERCES_CPP_NAMESPACE_QUALIFIER 
#endif
                         DOMElement* modeElt);

    /**
       Does work for public routines of this name
     */
    AcqSummaryInfo* getAcqSummaryWhere(const std::string& where);

    /**
       Utility does the work for getConstituentInfo and getConstituentByFswId
     */
    ConstitInfo* getConstituentWhere(const std::string& where);


    // Get info about each parameter file which was input source for the
    // specified entry in FSW_inputs.  Append to vector parms.
    bool getParmsFromInput(const std::string& fswInputKey, 
                           std::vector<ParmOffline>& parms);

    /// Returns key of file registered in Ancillary if there is one; else 0.
    /// @a ancClassKey is ancillary class key as string
    unsigned resolveAncAliasByKey(const std::string& alias, 
                                  const std::string& ancClassKey);

    /// Throw RdbException for error (e.g., no such voteKey, no ancillary
    /// found for one alias)
    bool resolveAncAliases(std::vector<std::string>& ancKeys,
                               const std::string& voteKeyStr);

    /// Given container vote file, return keys of vote files currently
    /// associated with precinct aliases in the container.
    bool resolveVoteAliases(std::vector<std::string>& voteKeys,
                               const std::string& ctnKeyStr);


    /// Update internal cache of parameter class names 
    /// for use of getParmsFromInput
    void updateParmClassCache(); 

    /**
       Return false if
          - voteKeyStr doesn't refer to vote in db
          - voteKey does exist, but @arg goodStatus was true and 
            vote did not have status "CREATED"
     */
    bool voteExists(const std::string& voteKeyStr, bool goodStatus=true);


    bool voteIsUpToDate(const std::string& voteKeyStr, 
                        std::vector<unsigned>* pk=0, bool* isCtn=0, 
                        bool force=false);


    // Following doesn't seem to be implemented anywhere
    //    bool containerIsUpToDate(const std::string& voteKeyStr);

    rdbModel::Rdb* m_rdb;
    MoodConnection* m_mood;
    bool            m_ownConnection; // if true, will need to delete 
    std::string     m_archive;  // translation of MOOT_ARCHIVE env. var.
    bool            m_dbg;
    std::vector<std::string> m_parmClasses; // cache them here

    // Cache filter assoc info; only re-retrieve if necessary
    static FilterAssoc* m_assoc;
  };
}

#endif
