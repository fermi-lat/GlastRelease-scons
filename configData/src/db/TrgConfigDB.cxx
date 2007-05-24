
#include "configData/db/TrgConfigDB.h"
#include "configData/gem/TrgConfigParser.h"
#include <cassert>

TrgConfigDB::TrgConfigDB(LatcDB* lc):m_lc(lc),m_key(0){}

bool TrgConfigDB::updateKey(const unsigned int key){
  int retval=false;
  if (key!=m_key){
    retval=true;
    std::string bcast=m_lc->getFilename("latc_DFT",key);
    TrgConfigParser p;
    int error=p.parse(this,bcast.c_str());
    assert(error==0);
    m_key=key;
  }
  return retval;
}
