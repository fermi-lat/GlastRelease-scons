
#include "configData/db/TrgConfigDB.h"
#include "configData/gem/TrgConfigParser.h"

TrgConfigDB::TrgConfigDB(LatcDB* lc):m_lc(lc),m_key(0){}

int TrgConfigDB::updateKey(unsigned key){
  int retval=0;
  if (key!=m_key){
    std::string bcast=m_lc->getFilename("latc_DFT",key);
    TrgConfigParser p;
    retval=p.parse(this,bcast.c_str());
    assert(retval==0);
    m_key=key;
  }
  return retval;
}
