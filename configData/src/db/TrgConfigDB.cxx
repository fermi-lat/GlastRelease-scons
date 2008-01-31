
#include "configData/db/TrgConfigDB.h"
#include "configData/gem/TrgConfigParser.h"
#include <cassert>

TrgConfigDB::TrgConfigDB(LatcDB* lc):m_lc(lc),m_key(0),m_allowMissing(false){}

bool TrgConfigDB::updateKey(const unsigned int key){
  int retval=false;
  bool newStyleConfig = false;
  if (key!=m_key){
    retval=true;
    std::string bcast=m_lc->getFilename("latc_DFT",key);
    if (bcast.length() == 0) // New style config 
      bcast = m_lc->getFilename("latc_GEM_TRG_GEM", key);
    newStyleConfig = true;
    std::cout<<bcast<<std::endl;
    TrgConfigParser p(m_allowMissing);
    int error=p.parse(this,bcast.c_str());
    assert(error==0);
    if (newStyleConfig) {
      bcast = m_lc->getFilename("latc_GEM_TRG_ROI", key);
      error = p.parse(this, bcast.c_str());
      assert(error==0);
    }
    m_key=key;
  }
  return retval;
}

void TrgConfigDB::allowMissing(bool allow){
  m_allowMissing=allow;
}
