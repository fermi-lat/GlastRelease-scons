
#include "configData/db/MootVoteParser.h"
#include "mootCore/FileDescrip.h"

#include "xmlBase/Dom.h"
#include "xmlBase/XmlParser.h"
#include <xercesc/dom/DOMElement.hpp>

#include <TString.h> 
#include <TRegexp.h> 

namespace configData {
  
  std::vector<std::string>           MootVoteParser::s_pNameVector;
  std::map<std::string,PrecinctType> MootVoteParser::s_pNameMap;

  void MootVoteParser::fillMaps() {
    if ( s_pNameVector.size() > 0 ) return;
    addToMaps(Unknown,"Unknown");
    addToMaps(TRG_GEM,"TRG_GEM"); 
    addToMaps(TRG_ROI,"TRG_ROI"); 
    addToMaps(GNL_Mode,"GNL_Mode"); 
    addToMaps(GNL_Timing,"GNL_Timing");
    addToMaps(ACD_Mode,"ACD_Mode"); 
    addToMaps(ACD_Timing,"ACD_Timing"); 
    addToMaps(ACD_Bias,"ACD_Bias"); 
    addToMaps(ACD_Hld,"ACD_Hld"); 
    addToMaps(ACD_PHA,"ACD_PHA");		      
    addToMaps(ACD_Veto,"ACD_Veto");		      
    addToMaps(CAL_Mode,"CAL_Mode"); 
    addToMaps(CAL_Timing,"CAL_Timing"); 
    addToMaps(CAL_FLE,"CAL_FLE"); 
    addToMaps(CAL_FHE,"CAL_FHE"); 
    addToMaps(CAL_LAC,"CAL_LAC"); 
    addToMaps(CAL_ULD,"CAL_ULD");
    addToMaps(TKR_Mode,"TKR_Mode"); 
    addToMaps(TKR_Timing,"TKR_Timing"); 
    addToMaps(TKR_Strips,"TKR_Strips"); 
    addToMaps(TKR_Thresh,"TKR_Thresh");
    addToMaps(FSW,"FSW");
    addToMaps(LPA,"LPA");
    addToMaps(ACD_LCI,"ACD_LCI"); 
    addToMaps(CAL_LCI,"CAL_LCI"); 
    addToMaps(TKR_LCI,"TKR_LCI");
    addToMaps(LCI_GLOBAL_ACD,"LCI_GLOBAL_ACD"); 
    addToMaps(LCI_GLOBAL_CAL,"LCI_GLOBAL_CAL"); 
    addToMaps(LCI_GLOBAL_TKR,"LCI_GLOBAL_TKR");    
  }

  
  void MootVoteParser::addToMaps(PrecinctType pt, const char* pName) {
    s_pNameVector.push_back(pName);
    s_pNameMap[pName] = pt;
  }

  const std::string& MootVoteParser::precinctName( PrecinctType pt ) {
    fillMaps();
    return pt > s_pNameVector.size() ? s_pNameVector[0] : s_pNameVector[pt];
  }

  PrecinctType MootVoteParser::precinctType( const std::string& pName ) {
    fillMaps();
    std::map<std::string,PrecinctType>::const_iterator itr = s_pNameMap.find(pName);
    return itr != s_pNameMap.end() ? itr->second : Unknown;
  }
  
  void MootVoteParser::printToText( std::ostream& os, MOOT::ConfigInfo& config, bool brief ) {
    if ( brief ) {
      os << config.getKey() << ' ' << config.getName() << std::endl;
      return;
    }
    os << "Key:       " << config.getKey() << std::endl;
    os << "Name:      " << config.getName() << std::endl;
    os << "Status:    " << config.getStatus() << std::endl;
    os << "Active:    " << config.getActive() << std::endl;
    os << "Time:      " << config.getCreationTime() << std::endl;
    os << std::endl;
  }
  

  void MootVoteParser::mergeSetsByUnion( std::set<unsigned>& master, const std::set<unsigned>& input ) {
    for ( std::set<unsigned>::const_iterator itr = input.begin(); itr != input.end(); itr++) {
      master.insert(*itr);
    }
  }
  
  void MootVoteParser::mergeSetsByISect( std::set<unsigned>& master, const std::set<unsigned>& input ) {
    for ( std::set<unsigned>::iterator itr = master.begin(); itr != master.end(); itr++) {
      if ( input.count(*itr) == 0 ) {
	master.erase(itr);
      }
    }
  } 
  
  bool MootVoteParser::regexpMatch(const std::string& exp, const std::string& str) {
    TRegexp re(exp.c_str());
    TString st(str.c_str());
    return st.Contains(re);
  }

  bool MootVoteParser::matchVoteTag( const DOMElement& elem, const VoteTag& voteTag ) {
    std::vector<DOMElement*> children;
    std::vector<DOMElement*> matched;
    xmlBase::Dom::getDescendantsByTagName(&elem,voteTag.m_tagName,children);
    for ( std::vector<DOMElement*>::const_iterator itr = children.begin(); itr != children.end(); itr++ ) {
      bool isMatch(true);
      for ( std::map<std::string,std::string>::const_iterator itrAttr = voteTag.m_tagAttr.begin(); 
	    itrAttr != voteTag.m_tagAttr.end(); itrAttr++ ) {
	std::string attrVal = xmlBase::Dom::getAttribute( *itr, itrAttr->first );
	if ( ! regexpMatch( itrAttr->second, attrVal ) ) {
	  isMatch = false;
	  break;
	}	
      }
      if ( isMatch ) matched.push_back(*itr);
    }

    if ( matched.size() == 0 ) return false;
    if ( voteTag.m_tagValue.size() == 0 ) return true;

    for ( std::vector<DOMElement*>::const_iterator itr2 = matched.begin(); itr2 != matched.end(); itr2++ ) {
      std::vector<DOMElement*> lexical;
      lexical.push_back(const_cast<DOMElement*>(*itr2));
      xmlBase::Dom::getDescendantsByTagName(*itr2,"singleton",lexical,false);
      xmlBase::Dom::getDescendantsByTagName(*itr2,"broadcast",lexical,false);
      for ( std::vector<DOMElement*>::const_iterator itr3 = lexical.begin(); itr3 != lexical.end(); itr3++ ) {
	std::string txtVal = xmlBase::Dom::getTextContent(*itr3);
	if ( regexpMatch( voteTag.m_tagValue, txtVal ) ) {
	  return true;
	}
      }
    }
    return false;
  }

  bool MootVoteParser::matchVoteContent( const std::string& fileName, const std::list<const VoteTag*>& votes, bool matchAll   ) {

    static XmlParser* parser(0);
    if ( parser == 0 ) {
      parser = new XmlParser(true);      
    }
    DOMDocument* doc(0);
    DOMElement* topElt(0);
    try {
      doc = parser->parse( fileName.c_str() );
      topElt = doc->getDocumentElement();
    } catch (ParseException ex) {
      std::cout << "caught exception with message " << std::endl;
      std::cout << ex.getMsg() << std::endl;
      return false;
    }

    for ( std::list<const VoteTag*>::const_iterator itr = votes.begin(); itr != votes.end(); itr++ ) {
      bool isMatch = matchVoteTag(*topElt,*(*itr));
      if ( isMatch ) {
	if ( !matchAll ) return true;
      } else {
	if ( matchAll ) return false;
      }	      
    }
    // only got here if all the tags came out the same way, which is the same as the matchAll flag
    return matchAll;
  }


  MootVoteParser::MootVoteParser()
    :m_mq(new MOOT::MootQuery){
  }

  MootVoteParser::~MootVoteParser(){
    delete m_mq;
  }

  
  /// function to print a set of configs
  bool MootVoteParser::printConfigs( std::ostream& os, const std::set<unsigned>& configs, bool brief ) {
    for ( std::set<unsigned>::const_iterator itr = configs.begin(); itr != configs.end(); itr++ ) {
      MOOT::ConfigInfo* ci = m_mq->getConfigInfo(*itr);
      if ( ci == 0 ) return false;
      printToText(os,*ci,brief);
    }
    return true;
  }

  /// function to get all the configs built from a set of votes
  bool MootVoteParser::getConfigsFromVotes( std::set<unsigned>& configs, const std::set<unsigned>& votes ) {    
    for ( std::set<unsigned>::const_iterator itr = votes.begin(); itr != votes.end(); itr++ ) {
      if ( ! m_mq->getConfigsForVote(*itr,configs) ) return false;
    }
    return true;
  }
  

  /// function to retrieve filename for given key key
  bool MootVoteParser::getVoteFilename(std::string& fileName, unsigned key) {
    MOOT::VoteInfo* vi = m_mq->getVoteInfo(key);
    if ( vi == 0 ) return false;
    fileName = std::string(getenv("MOOT_ARCHIVE"));
    fileName += "/";
    fileName += vi->getSrc();
    return true;
  }

  /// function to get all the master keys that match precinct vote alias
  bool MootVoteParser::getMasterKeysWithVoteAliases( std::set<unsigned>& keys, const std::list<VoteAlias>& votes, 
						     bool matchAll ) {
    bool first(true);
    for ( std::list<VoteAlias>::const_iterator itr = votes.begin(); itr != votes.end(); itr++ ) {
      std::set<unsigned> tempSet;
      const std::string pn = precinctName(itr->m_precinct);
      if ( ! m_mq->getConfigsForAlias(pn,itr->m_alias,tempSet) ) return false;
      if ( first ) {
	first = false;
	keys = tempSet;
      } else {
	if ( matchAll ) {
	  mergeSetsByUnion(keys,tempSet);
	} else {
	  mergeSetsByISect(keys,tempSet);
	}
      }
    }
    return true;
  }

  /// function to get all the master keys that match precinct vote content
  bool MootVoteParser::getMasterKeysWithVoteContent( std::set<unsigned>& keys, const std::list<VoteTag>& voteTags, 
						     bool matchAll ) {

    std::map<PrecinctType,std::vector< std::pair<unsigned, std::string> > > voteIdMap;
    std::map<PrecinctType,std::list<const VoteTag*> > voteTagMap;
    for ( std::list<VoteTag>::const_iterator itr = voteTags.begin(); itr != voteTags.end(); itr++ ) {
      voteTagMap[itr->m_precinct].push_back(&(*itr));
      std::vector<std::pair<unsigned, std::string> >& votes = voteIdMap[itr->m_precinct];      
      if ( votes.size() == 0 ) {	
	std::vector<unsigned> voteKeys;
        m_mq->listVoteKeys( voteKeys, precinctName(itr->m_precinct) );
	for ( std::vector<unsigned>::const_iterator itrVk = voteKeys.begin(); itrVk != voteKeys.end(); itrVk++ ) {
	  std::string voteFileName;
	  if ( ! getVoteFilename(voteFileName,*itrVk) ) return false;
	  votes.push_back( std::make_pair( *itrVk, voteFileName) );
	}
      }
    }

    bool first(true);
    for ( std::map<PrecinctType,std::list<const VoteTag*> >::const_iterator itr2 = voteTagMap.begin();
	  itr2 != voteTagMap.end(); itr2++ ) {
      std::set<unsigned> voteSet;
      std::set<unsigned> tempSet;
      const std::vector< std::pair<unsigned, std::string> >&  votes = voteIdMap[itr2->first];
      for ( std::vector<std::pair<unsigned, std::string> >::const_iterator itrV = votes.begin(); itrV != votes.end(); itrV++ ) {
	bool isMatch = matchVoteContent(itrV->second,itr2->second,matchAll);
	if ( isMatch ) {	  
	  voteSet.insert(itrV->first);
	}
      }
      if ( ! getConfigsFromVotes(tempSet,voteSet ) ) return false;
      if ( first ) {
	first = false;
	keys = tempSet;
      } else {
	if ( matchAll ) {
	  mergeSetsByUnion(keys,tempSet);
	} else {
	  mergeSetsByISect(keys,tempSet);
	}
      }
    } 
    return true;
  }
  
  /// function to get all the master keys that match ancillary files
  bool MootVoteParser::getMasterKeysWithAncilaryFiles( std::set<unsigned>& keys, const std::list<unsigned>& ancilKeys, 
						       bool matchAll ) {
    bool first(true);
    for ( std::list<unsigned>::const_iterator itr = ancilKeys.begin(); itr != ancilKeys.end(); itr++ ) {
      std::set<unsigned> tempSet;
      if ( ! m_mq->getConfigsForAncillary(*itr, tempSet) ) return false;
      if ( first ) {
	first = false;
	keys = tempSet;
      } else {
	if ( matchAll ) {
	  mergeSetsByUnion(keys,tempSet);
	} else {
	  mergeSetsByISect(keys,tempSet);
	}
      }
    }
    return true;
    
  }   


}
