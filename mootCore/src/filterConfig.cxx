//  $Header$
/**
       @file filterConfig.cxx
       That part of implementation of MootQuery class concerning filter config.
       Includes MootQuery methods and implementation of helper class
       FilterAssoc
*/
// stl
#include <string>
#include <cstdio>
#include <stdexcept>
#include <vector>
#include <algorithm>
#include <map>   // I think we'll need it; TBD

// Our includes
#include "mootCore/MootQuery.h"
#include "mootCore/MoodConnection.h"
#include "mootCore/DbUtil.h"
#include "mootCore/FileDescrip.h"

#include "xmlBase/XmlParser.h"
#include "xmlBase/XmlErrorHandler.h"
#include "xmlBase/Dom.h"

#include "facilities/Util.h"

#ifndef SWIG
XERCES_CPP_NAMESPACE_BEGIN
#endif
class  DOMDocument;
class  DOMElement;
#ifndef SWIG
XERCES_CPP_NAMESPACE_END
#endif

// external library include
// from external libraries
#include <xercesc/dom/DOMElement.hpp>
#include <xercesc/dom/DOMDocument.hpp>

namespace MOOT {
  XERCES_CPP_NAMESPACE_USE

#define INVALID_MODE 0xffffffff

  //  typedef std::pair<ConstitInfo*, unsigned> FilterInfo; // 
  struct FilterInfo {
    FilterInfo(ConstitInfo* pI=0, unsigned id=0xffffffff, 
               const std::string& name=std::string("")) :
      pInfo(pI), hId(id), hName(name) {}
    // copy constructor
    FilterInfo(const FilterInfo& f) : pInfo(f.pInfo),
                                      hId(f.hId), hName(f.hName) {}
    ~FilterInfo() {}
    ConstitInfo* pInfo;
    unsigned     hId;
    std::string  hName;
  };

  typedef std::vector<FilterInfo> AcqMode;
  class FilterAssoc {
  public:
    FilterAssoc() : m_configKey(0)  {initModeNames();}

    ~FilterAssoc() {
      clean(); 
      m_modeNames->clear();
    }

    /// Add new constituent info to store; return pointer to it
    ConstitInfo* addConstit(ConstitInfo* pConstit);

    void clean();

    /// Search m_constits for entry with constituent name = @a name
    ConstitInfo* findConstit(const std::string& name);


    /// Find index in static array of @a modeName
    unsigned modeIndex(const std::string& modeName) const;

    std::vector<ConstitInfo*> m_constits;

    unsigned m_configKey;
    std::string m_parmKeyStr;
    std::string m_sbsKeyStr;        // fmx export key
    std::string m_sbsFkStr;         // MOOT key in FSW_inputs table
    std::string m_lpaDbName;

    AcqMode  m_modeData[LPA_MODE_count];
  private:
    void initModeNames();

    static std::vector<std::string>* m_modeNames;
  };

  // Initialize a couple statics
  FilterAssoc* MootQuery::m_assoc = 0;
  std::vector<std::string>*  FilterAssoc::m_modeNames = 0;


  // Some MootQuery methods
  bool MootQuery::cacheFilterConfig(unsigned configKey) {
    using facilities::Util;

    if (!configKey) return false;
    if (!m_assoc) m_assoc = new FilterAssoc();
    else if (configKey == m_assoc->m_configKey) return true;

    // Fetch key, (relative) src path for parameter file of type 
    std::string configKeyStr;
    Util::utoa(configKey, configKeyStr);
    std::string parmClassKeyStr = 
      DbUtil::getColumnValue(m_rdb, "Parameter_class", "Parameter_class_key",
                             "name", "fsw_FSW_Filter_Associate");

    if (parmClassKeyStr.empty() ) {
      return false;  // and maybe do other stuff
    }

    std::string where(" WHERE config_fk='");
    where += configKeyStr + 
      std::string("' and Parameter_fk=Parameters.parm_key and ") +
      std::string("Parameters.class_fk='") + parmClassKeyStr +
      std::string("'");
    std::string parmKeyStr =  
      DbUtil::getColumnWhere(m_rdb, "Configs_to_Parameters,Parameters",
                             "parm_key", where, false);

    // check to make sure we got something
    if (parmKeyStr.empty() ) {
      return false;  // and maybe do other stuff
    }

    // Similarly for sbs key
    std::string sbsKeyStr = DbUtil::getColumnValue(m_rdb, "Configs",
                                                   "sbs_key", "config_key",
                                                   configKeyStr);
    if (sbsKeyStr.empty() ) return false;

    std::string sbsFkStr = DbUtil::getColumnValue(m_rdb, "FSW_inputs",
                                                  "FSW_input_key", 
                                                  "FSW_id", sbsKeyStr);

    if ((parmKeyStr == m_assoc->m_parmKeyStr) && 
        (sbsKeyStr == m_assoc->m_sbsKeyStr) )  {
      m_assoc->m_configKey = configKey;
      return true;
    }

    // Old information is no good; 
    if (m_assoc) m_assoc->clean();
    else m_assoc = new FilterAssoc();
    m_assoc->m_configKey = configKey;
    m_assoc->m_parmKeyStr = parmKeyStr;
    m_assoc->m_sbsKeyStr = sbsKeyStr;
    m_assoc->m_sbsFkStr = sbsFkStr;
  

    // Make abs. path
    std::string path = m_archive + 
      DbUtil::getColumnValue(m_rdb, "Parameters", "source", "parm_key",
                             parmKeyStr);
    if (path == m_archive) {
      return false;
    }

    // Parse
    return parseFilterParm(path);
        
  }

  // Public routines
  unsigned MootQuery::getActiveFilters(unsigned configKey, 
                                       std::vector<ConstitInfo>& info,
                                       LpaMode mode, bool clear) {
    if (clear) info.clear();
    if (!cacheFilterConfig(configKey)) return 0;
    unsigned cnt;

    if (mode == LPA_MODE_ALL) {
      cnt = m_assoc->m_constits.size();
      for (unsigned ix = 0; ix < cnt; ix++) {
        info.push_back(*((m_assoc->m_constits)[ix]) );
      }
      return cnt;
    }
    else if (mode >= LPA_MODE_count) return 0;
    
    // legitimate mode
    AcqMode& modeData = m_assoc->m_modeData[mode];
    cnt = modeData.size();            
    for (unsigned fIx = 0; fIx < cnt; fIx++) {
      info.push_back(*((modeData[fIx]).pInfo));
    }
    return cnt;
  }
      
  const ConstitInfo* MootQuery::getActiveFilter(unsigned configKey, 
                                                LpaMode mode,
                                                unsigned handlerId,
                                                std::string& handlerName) {
    handlerName.clear();
    if (!cacheFilterConfig(configKey)) return 0;
    if (mode >= LPA_MODE_count) return 0;   // unrecognized mode
    // attempt to find entry matchin handlerId for this mode
    AcqMode& modeData = m_assoc->m_modeData[mode];
    unsigned cnt = modeData.size();            

    for (unsigned fIx = 0; fIx < cnt; fIx++) {
      if (modeData[fIx].hId == handlerId) {
        handlerName = modeData[fIx].hName;
        return modeData[fIx].pInfo;
      }
    }
    return 0;
  }



  bool MootQuery::parseFilterParm(const std::string& path) {
    using xmlBase::Dom;

    xmlBase::XmlParser  parser;
    DOMDocument* doc = parser.parse(path.c_str());
    if ( doc == 0 )  return 0;

    DOMElement* elt = doc->getDocumentElement();

    if (!Dom::checkTagName(elt, "Associations")) {
      std::cerr << path << " is not a valid Filter_Associate vote" 
                << std::endl;
      return false;
    }
    elt = Dom::findFirstChildByName(elt, "LPA_DB");
    if (!elt) {
      std::cerr << path << " is not a valid Filter_Associate vote" 
                << std::endl;
      return false;
    }
    m_assoc->m_lpaDbName = Dom::getAttribute(elt, "name");
    std::vector<DOMElement*> modeElts;

    Dom::getChildrenByTagName(elt, "LPA_Mode", modeElts);
    for (unsigned ix = 0; ix < modeElts.size(); ix++) {
      if (!parseLpaModeElt(modeElts[ix])) return false;
    }
    return true;
  }


  bool MootQuery::parseLpaModeElt(DOMElement* modeElt) {
    using xmlBase::Dom;

    std::string modeString = Dom::getAttribute(modeElt, "name");
    if (modeString.empty()) return true; // nothing to do

    std::vector<std::string> tokens;
    facilities::Util::stringTokenize(modeString, " ", tokens);
    unsigned firstIx = m_assoc->modeIndex(tokens[0]);

    if (firstIx == INVALID_MODE) {
      std::cerr << "Invalid Filter_Associate parameter. Unknown mode name " 
                << tokens[0] << std::endl;
      return false;
    }
    DOMElement* activeElt = Dom::findFirstChildByName(modeElt, "Filters");
    if (!activeElt) return false;

    std::string activeList = Dom::getAttribute(activeElt, "active");
    std::vector<std::string> activeTokens;

    facilities::Util::stringTokenize(activeList, " ", activeTokens);

    AcqMode& firstMode=m_assoc->m_modeData[firstIx];

    if (!firstMode.empty() ) {
      std::cerr << "Invalid Filter_Associate parameter. Duplicate mode"
                << tokens[0] << std::endl;
      return false;
    }
    std::vector<DOMElement*> handlerElts;
    Dom::getChildrenByTagName(modeElt, "Handler", handlerElts);
    

    firstMode.reserve(handlerElts.size());
    for (unsigned hIx = 0; hIx < handlerElts.size(); hIx++) {
      std::string handName = Dom::getAttribute(handlerElts[hIx], "name");
      
      if (find(activeTokens.begin(), activeTokens.end(), handName) == 
          activeTokens.end() )  continue;

      std::string instance = Dom::getAttribute(handlerElts[hIx], "instance");
      std::string handIdStr = Dom::getAttribute(handlerElts[hIx], "id");
      if (instance.empty() || handIdStr.empty()) {
        std::cerr << "Invalid Filter_Associate parameter; null reqd attribute"
                  << std::endl;
        return false;
      }
      unsigned handId = facilities::Util::stringToUnsigned(handIdStr);
      ConstitInfo* pInfo = m_assoc->findConstit(instance);
      if (!pInfo) {
        std::string where(" WHERE sbs_fk='");
        where += m_assoc->m_sbsFkStr + 
          std::string("' AND Constituents.prim_key=Constituents_fk AND ");
        where += std::string("Constituents.name='") + instance + 
          std::string("'");
        std::string constitKeyStr = 
          DbUtil::getColumnWhere(m_rdb, "SBS_to_Constituents,Constituents",
                                 "Constituents.prim_key", where, false);
        if (constitKeyStr.empty()) return false;
        
        where = std::string(" WHERE prim_key='" );
        where += constitKeyStr + std::string("'");
        pInfo = getConstituentWhere(where);
        if (!pInfo) return false;
        pInfo = m_assoc->addConstit(pInfo);
        //        delete pNewInfo;
      }
      FilterInfo f(pInfo, handId, handName);
      firstMode.push_back(f);
    }   // finishes processing of first mode in the list.  
    // For remaining modes, just copy
    for (unsigned tokIx = 1; tokIx < tokens.size(); tokIx++) {
      unsigned modeIx = m_assoc->modeIndex(tokens[tokIx]);
      if (modeIx == INVALID_MODE) {
        std::cerr << "Invalid Filter_Associate parameter. Unknown mode name " 
                  << tokens[tokIx] << std::endl;
        return false;
      }
      AcqMode& curMode = m_assoc->m_modeData[modeIx];
      if (!curMode.empty() ) {
        std::cerr << "Invalid Filter_Associate parameter. Duplicate mode"
                  << tokens[tokIx] <<  std::endl;
        return false;
      }
      curMode.reserve(handlerElts.size());
      curMode = firstMode;
    }
    return true;
  }


  // FilterAssoc implementation
  ConstitInfo* FilterAssoc::addConstit(ConstitInfo* pConstit) {
    m_constits.push_back(pConstit);
    return pConstit;
  }

  void FilterAssoc::clean() {
    for (unsigned ix=0; ix < m_constits.size(); ix++) {
      delete m_constits[ix];
    }
    m_constits.clear();
    for (unsigned ix = 0; ix < LPA_MODE_count; ix++) {
      m_modeData[ix].clear();
    }
    m_configKey = 0;
    m_parmKeyStr.clear();
    m_lpaDbName.clear();
  }

  ConstitInfo* FilterAssoc::findConstit(const std::string& name)  {
    for (unsigned ix = 0; ix < m_constits.size(); ix++) {
      if (name == (m_constits[ix])->getName()) return m_constits[ix];
    }
    return 0;
  }


  #define N_MODES 8
  void FilterAssoc::initModeNames() {
    if (!m_modeNames) {
      m_modeNames = new std::vector<std::string>();
      m_modeNames->reserve(N_MODES);
    }
    if (m_modeNames->size() == 0) {
      m_modeNames->push_back("NORM");
      m_modeNames->push_back("TOO");
      m_modeNames->push_back("GRB0");
      m_modeNames->push_back("GRB1");
      m_modeNames->push_back("GRB2");
      m_modeNames->push_back("SOLAR");
      m_modeNames->push_back("CALIB");
      m_modeNames->push_back("DIAG");
    }
  }

  unsigned FilterAssoc::modeIndex(const std::string& modeName) const {
    for (unsigned ix = 0; ix < m_modeNames->size(); ix++) {
      if (modeName == (*m_modeNames)[ix]) {
        return ix;
      }
    }
    return INVALID_MODE;
  }

}
