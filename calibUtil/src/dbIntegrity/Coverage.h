// $Header$
#ifndef CALIBUTIL_COVERAGE_H
#define CALIBUTIL_COVERAGE_H
/**
    @file Coverage.h
*/
#include <vector>
#include <string>
#include "facilities/Timestamp.h"

namespace rdbModel {
  class Rdb;
  class Connection;
}
namespace calibUtil {
  class Metadata;
}

/** 
    @class Coverage
    
    Check whether calibrations for a particular calibration type (or each in
    a predefined collection) exist covering full time interval, and whether
    they overlap.
*/
class Coverage {
public:
  Coverage(calibUtil::Metadata* meta,
           const std::string& instr, const std::string& flavor, 
           const std::string& level, const facilities::Timestamp& ts);

  /// The @nickname argument may be a single legal calib_type or may be
  /// one of the special wildcard values "CAL" "TKR" "*", in which case it
  /// is expanded to a suitable list of calib types (if we have a schema) 
  bool expandTypes(std::string& nickname, std::vector<std::string>& types);

  ~Coverage() {};
  
  // raison d'etre of this class.  Examine all db rows of type @calibType
  // and also matching flavor, instrument, level, etc., to see if there
  // is precisely one valid an any one time (with some tolerance for
  // small overlap)
  unsigned checkType(std::string calibtype);
  void setOverlap(unsigned seconds) {m_overlap = seconds;}
private:
  rdbModel::Connection* m_conn;
  rdbModel::Rdb* m_rdb;
  std::string m_instr;
  std::string m_flavor;
  std::string m_level;
  std::string m_calib;

  facilities::Timestamp m_ts;
  std::string    m_table;
  int m_overlap;
  std::vector<std::string> m_selects;
  std::vector<std::string> m_orderBy;

};
#endif
