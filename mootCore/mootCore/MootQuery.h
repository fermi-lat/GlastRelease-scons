// $Header$
// Handles registering a single file into Parameters table.  
// Might also handle other parts of building a config; TBD.
#ifndef MOOT_MootQuery_h
#define MOOT_MootQuery_h

#include <string>
#include <vector>
#include "rdbModel/Tables/Column.h"
#include "mootCore/FileDescrip.h"

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

typedef std::vector<std::string> VectorOfString;
namespace MOOT {
  class MoodConnection;

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

    /** Store everything of interest about this row in AncillaryAlias
        table in AncAliasInfo structure.
        @return 0 if no row exists with this key;  otherwise ptr
         to filled structure.  Caller is responsible for deleting.
     */
    AncAliasInfo* getAncAliasInfo(unsigned key);

    unsigned getAncClasses(std::vector<std::string>& names);

    /** Store everything of interest about this row in Ancillary
        table in AncInfo structure.
        @return 0 if no row exists with this key; otherwise ptr
         to filled structure.  Caller is responsible for deleting.
     */
    AncInfo* getAncInfo(unsigned key);

    // **TODO**


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
       alg_step = step.  Empty string for @status or @instr argument
       means any value is acceptable.  Similarly, @step=0 means
       Config can have any value for step.
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
       Return keys of parameter entries supplied when config was
       created.
     */
    bool getConfigParmsRequest(unsigned configKey, 
                               std::vector<unsigned>& parameterKeys);
    /**
       Return keys of parameter entries used to build FSW inputs
       for this config
     */
    bool getConfigParmsUsed(unsigned configKey, 
                            std::vector<unsigned>& parameterKeys);


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

    /** Store everything of interest about this row in Parameters
        table in Info structure.
        @return 0 if no row exists with this key; 1 otherwise
     */
    ParmInfo* getParmInfo(unsigned key);


    /// Fill supplied argument with parameter class names. Return count
    unsigned getParmClasses(std::vector<std::string>& names);

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

    /* If vote file specified by @a voteKey is up to date, put keys 
       of parameter files of all relevant classes associated with the
       vote file in the @a parmKeys vector and return true.  
       Else return false and an empty list.
     */
    bool getVoteParameters(unsigned voteKey, 
                           std::vector<unsigned>& parmKeys);

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

    /// The "vote" referred to may concern a single precinct or may
    /// be container file, referring to several individual precinct files.
    /// If optional arg. is supplied and vote is up to date, return
    /// parm keys associated with the vote
    bool voteIsUpToDate(unsigned voteKey, std::vector<unsigned>* parmKeys=0);

  private:
    /*
    unsigned getVoteParmViaClassKey(const std::string& voteKeyStr,
                                    const std::string& parameterClassKey);
    */


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

    bool voteIsUpToDate(const std::string& voteKeyStr, 
                        std::vector<unsigned>* pk=0);

    bool containerIsUpToDate(const std::string& voteKeyStr);

    rdbModel::Rdb* m_rdb;
    MoodConnection* m_mood;
    bool            m_ownConnection; // if true, will need to delete 
    std::string     m_archive;  // translation of MOOT_ARCHIVE env. var.
    bool            m_dbg;
    std::vector<std::string> m_parmClasses; // cache them here
  };
}

#endif
