
#include "configData/db/TrgConfigDB.h"
#include "configData/gem/TrgConfigParser.h"
#include <cassert>

TrgConfigDB::TrgConfigDB(LatcDB* lc):m_lc(lc),m_key(0),m_allowMissing(false){}

bool TrgConfigDB::updateKey(const unsigned int key){
  int retval=false;
  if (key!=m_key){
    TrgConfigParser parser;
    retval=true;
    std::string bcast=m_lc->getFilename("latc_DFT",key);
    if (bcast.length() != 0){ // Old style config 
      int error=parser.parse(this,bcast.c_str());
      std::cout<<bcast<<std::endl;
      assert(error==0);
    } else { // new style configuration
      bcast = m_lc->getFilename("latc_DFT_TRG_GEM", key); //this one may or may not exist
      if (bcast.length()!=0){
        int error = parser.parse(this, bcast.c_str());
        std::cout<<bcast<<std::endl;
        assert(error==0);
      }
      bcast = m_lc->getFilename("latc_GEM_TRG_GEM", key);
      int error=parser.parse(this,bcast.c_str());
      std::cout<<bcast<<std::endl;
      assert(error==0);
      bcast = m_lc->getFilename("latc_GEM_TRG_ROI", key);
      error = parser.parse(this, bcast.c_str());
      std::cout<<bcast<<std::endl;
      assert(error==0);
    }
    if(!m_allowMissing&& !parser.configComplete()){
      std::cerr<<"The configuration is not complete. Missing tags: ";
      std::list<std::string> tags=parser.remainingTags();
      for (std::list<std::string>::iterator it=tags.begin();it!=tags.end();it++)std::cerr<<(*it)<<" ";
      std::cerr<<std::endl;
      assert(parser.configComplete());
    }
    m_key=key;
  }
  return retval;
}

void TrgConfigDB::allowMissing(bool allow){
  m_allowMissing=allow;
}
