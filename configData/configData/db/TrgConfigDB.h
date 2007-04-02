/** @class TrgConfigDB
    @author M. Kocian

    A simple DB class for TrgConfig

    $Id$
*/
#ifndef TRGCONFIGDB_H
#define TRGCONFIGDB_H

#include "configData/gem/TrgConfig.h"
#include "configData/db/LatcDB.h"

class TrgConfigDB:public TrgConfig{
 public:
  /// Constructor with the DB interface implementation as the argument
  TrgConfigDB(LatcDB* lc);
  /// Check if we need to parse a new object, parse it if necessary
  int updateKey(unsigned key);
 private:
  LatcDB* m_lc;
  unsigned m_key;
};
    
  
#endif
