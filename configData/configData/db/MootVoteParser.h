//      -*- Mode: c++ -*-
/** @class MootVoteParser
    @author E. Charles
    

    Class to get configuratins from MOOT that match certain Intents

    $Id$
*/
#ifndef MOOTVOTEPARSER_H
#define MOOTVOTEPARSER_H

#include <string>
#include <vector>
#include <set>
#include <map>
#include <list>
#include <iostream>
#include "mootCore/MootQuery.h"
#include "mootCore/FileDescrip.h"

// XML
#include "xmlBase/Dom.h"
#include <xercesc/dom/DOMElement.hpp>


namespace configData {

  using namespace xmlBase;

  enum PrecinctType { Unknown,
                      TRG_GEM, TRG_ROI, 
		      GNL_Mode, GNL_Timing,
		      ACD_Mode, ACD_Timing, ACD_Bias, ACD_Hld, ACD_PHA, ACD_Veto,
		      CAL_Mode, CAL_Timing, CAL_FLE, CAL_FHE, CAL_LAC, CAL_ULD,
		      TKR_Mode, TKR_Timing, TKR_Strips, TKR_Thresh,
		      FSW,
		      LPA,
		      ACD_LCI, CAL_LCI, TKR_LCI, 
		      LCI_GLOBAL_ACD, LCI_GLOBAL_CAL, LCI_GLOBAL_TKR };


  struct VoteAlias {
  public:
    VoteAlias(PrecinctType pt, const std::string& alias)
      :m_precinct(pt),m_alias(alias){;}
    VoteAlias():m_precinct(Unknown){;}
    ~VoteAlias(){;}
    VoteAlias(const VoteAlias& other)
      :m_precinct(other.m_precinct),
       m_alias(other.m_alias){;}
    PrecinctType   m_precinct;
    std::string    m_alias;
  };
  
  struct VoteTag {
  public:
    VoteTag(PrecinctType pt, const std::string& tagName, const std::string tagValue)
      :m_precinct(pt),m_tagName(tagName),m_tagValue(tagValue){;}
    VoteTag():m_precinct(Unknown){;}
    ~VoteTag(){;}
    VoteTag(const VoteTag& other)
      :m_precinct(other.m_precinct),
       m_tagName(other.m_tagName),
       m_tagValue(other.m_tagValue),
       m_tagAttr(other.m_tagAttr){;}
    PrecinctType                        m_precinct;
    std::string                         m_tagName;
    std::string                         m_tagValue;
    std::map<std::string,std::string>   m_tagAttr;       
  };


  class MootVoteParser {
    
  public:
 
    static const std::string& precinctName( PrecinctType pt );    

    static PrecinctType precinctType( const std::string& pName );    

    static void printToText( std::ostream& os, MOOT::ConfigInfo& config, bool brief  );

    static void mergeSetsByUnion( std::set<unsigned>& master, const std::set<unsigned>& input );

    static void mergeSetsByISect( std::set<unsigned>& master, const std::set<unsigned>& input );

    static bool regexpMatch(const std::string& exp, const std::string& str);

    static bool matchVoteTag(const DOMElement& elem, const VoteTag& voteTag);

    static bool matchVoteContent( const std::string& fileName, const std::list<const VoteTag*>& vote, bool matchAll  );

  private:

    static void fillMaps();
    static void addToMaps(PrecinctType pt, const char* pName);

    static std::vector<std::string>           s_pNameVector;
    static std::map<std::string,PrecinctType> s_pNameMap;

  public:
    
    /// Default constructor
    MootVoteParser();

    /// D'tor
    ~MootVoteParser();

    /// function to print a set of configs
    bool printConfigs( std::ostream& os, const std::set<unsigned>& configs, bool brief );

    /// function to get all the configs built from a set of votes
    bool getConfigsFromVotes( std::set<unsigned>& configs, const std::set<unsigned>& votes );

    /// function to retrieve filename for given key key.
    bool getVoteFilename(std::string& fileName, unsigned key);
    
    /// function to get all the master keys that match precinct vote alias
    bool getMasterKeysWithVoteAliases( std::set<unsigned>& keys, const std::list<VoteAlias>& votes, bool matchAll );

    /// function to get all the master keys that match precinct vote content
    bool getMasterKeysWithVoteContent( std::set<unsigned>& keys, const std::list<VoteTag>& voteTags, bool matchAll );

    /// function to get all the master keys that match ancillary files
    bool getMasterKeysWithAncilaryFiles( std::set<unsigned>& keys, const std::list<unsigned>& ancilKeys, bool matchAll );      

  private:
    
    MOOT::MootQuery* m_mq;
    
  };

};

#endif
