
#include "configData/db/TrgConfigDB.h"
#include "configData/gem/TrgConfigParser.h"
#include <cassert>

TrgConfigDB::TrgConfigDB(LatcDB* lc):m_lc(lc),m_key(0),m_allowMissing(false){}

bool TrgConfigDB::updateKey(const unsigned int key){
  int retval=false;
  if (key!=m_key){
    TrgConfigParser p(m_allowMissing);
    retval=true;
    std::string bcast=m_lc->getFilename("latc_DFT",key);
    if (bcast.length() != 0){ // Old style config 
      int error=p.parse(this,bcast.c_str());
      std::cout<<bcast<<std::endl;
      assert(error==0);
    } else { // new style configuration
      bcast = m_lc->getFilename("latc_DFT_TRG_GEM", key); //this one may or may not exist
      if (bcast.length()!=0){
        int error = p.parse(this, bcast.c_str());
        std::cout<<bcast<<std::endl;
        assert(error==0);
      }
      bcast = m_lc->getFilename("latc_GEM_TRG_GEM", key);
      int error=p.parse(this,bcast.c_str());
      std::cout<<bcast<<std::endl;
      assert(error==0);
      bcast = m_lc->getFilename("latc_GEM_TRG_ROI", key);
      error = p.parse(this, bcast.c_str());
      std::cout<<bcast<<std::endl;
      assert(error==0);
    }
    m_key=key;
  }
  return retval;
}

void TrgConfigDB::allowMissing(bool allow){
  m_allowMissing=allow;
}
