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
  /// Virtual destructor
  virtual ~TrgConfigDB(){}
  /// Check if we need to parse a new object, parse it if necessary
  bool updateKey(const unsigned int key);
  /// Allow partial configurations to be used
  void allowMissing(bool allow);
  /// Allowing partial configs?
  bool allowingMissing(){return m_allowMissing;}
  
 protected:
  LatcDB* m_lc;
  unsigned m_key;
  bool m_allowMissing;
};
    
  
#endif
