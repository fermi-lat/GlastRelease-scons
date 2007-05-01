// $Header$
// Handles registering a single file into Parameters table.  
// Might also handle other parts of building a config; TBD.
#ifndef MOOT_DbUtil_h
#define MOOT_DbUtil_h

#include <string>

/** @file MootDbUtil.h
    @author J. Bogart
  
    Declare the class MootDbUtil. Includes services for
    common low-level query operations.  This class is
    only for use within the MOOT package.

*/
namespace rdbModel {
  class Rdb;
  class ResultHandle;
}


namespace MOOT {
  class DbUtil   {
  public:
    static int checkResults(rdbModel::ResultHandle* res,
                            const std::string& errMsgPrefix, 
                            int expected);

    /// Utility to get a single column value, given table name and key
    /// Optionally declare intent to update
    static std::string getColumnValue(rdbModel::Rdb* rdb, 
                                      const std::string& table, 
                                      const std::string& col,
                                      const std::string& keyCol,
                                      const std::string& key,
                                      bool forUpdate=false);
    /// Utility to get a single column value, given table name and 
    /// valid 'where' clause (should start " WHERE.." )
    /// If 'onlyOne' is true (default), complain if more - or less-
    ///  are found
    /// Else just return first one.
    static std::string getColumnWhere(rdbModel::Rdb* rdb, 
                                      const std::string& table, 
                                      const std::string& col,
                                      const std::string& where, 
                                      bool onlyOne=true,
                                      bool forUpdate=false);

    /// Utility to get a single column value (ascending order) 
    /// from possibly multiple
    /// rows.   Return value is number of rows found, or -1 for error.
    /// This function will append results to vals vector
    static int getAllWhere(rdbModel::Rdb* rdb,
                            const std::string& table,
                            const std::string& col,
                            const std::string& where,
                            std::vector<std::string>& vals);

    /**
       return keys (up to value of limit if limit > 0)  for all 
       rows satisfying condition in where.  By default return in descending
       order, so newest is always returned. May specify ascending, 
       however.
     */
    static unsigned getKeys(std::vector<unsigned>& keys,
                            rdbModel::Rdb* rdb, const std::string& table,
                            const std::string& keyCol, 
                            const std::string& where, unsigned limit=0,
                            bool ascend=false);
                            

    /** Utility to update value of a single column for rows satisfying
        the condition in 'where'.  Returns number of rows updated, or
        -1 for failure.
     */
    static int setColumnValue(rdbModel::Rdb* rdb,
                              const std::string& table,
                              const std::string& colName,
                              const std::string& value,
                              const std::string& where);
 /*, bool onlyOne=true); What was this for?? */

    /** Utility to update value of a single column for row satisfying
        the condition of keyCol = key.  Returns number of rows updated
        (should be 1, or possibly 0) or  -1 for failure.
     */
    static int setColumnValue(rdbModel::Rdb* rdb,
                              const std::string& table,
                              const std::string& colName,
                              const std::string& value,
                              const std::string& keyCol,
                              const std::string& key);

    static bool startTransaction(rdbModel::Rdb* rdb);
    static bool commit(rdbModel::Rdb* rdb);    

  private:

  };
}

#endif
