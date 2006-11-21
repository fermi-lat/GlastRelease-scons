// $Header$
// Thin wrapper for rdb connection
#ifndef MOOT_MoodConnection_h
#define MOOT_Connection_h

#include <string>
#include <vector>

/** @file MoodConnection.h
    @author J. Bogart
  
    Declare the class MoodConnection which maintains connection to MOOD
    database.
    
    For write access, must have a  ~/.my.cnf with user, password for
    group MOOT_WRITE.  A typical file might contain lines like this:

        [MOOT_WRITE]
        user=user_with_priv
        password="write_pwd"

    where some_userid, the_password, user_with_priv and write_pwd are 
    replaced with legitimate values.

    ~/.my.cnf should only be readable by the owner of the file.

*/
namespace rdbModel {
  class Rdb;
  class Connection;
}


namespace MOOT {

  class MoodConnection   {
  public:
    /// @arg wrt, input.  Does client need write access to db?
    /// There is an implicit argument:  environment variable MOOT_ARCH
    /// must have a sensible setting
    MoodConnection(bool wrt=false, bool dbg=false);

    rdbModel::Rdb* getRdb() {return m_rdb;}

    bool isWrite() const {return m_write;}

    bool isDebug() const {return m_dbg;}
    std::string getNickname() const {return m_nick;}

    ~MoodConnection();
    
  
  private:
    
    // Open connection to db, check that db schema (embodied in xml
    // file) is compatible.
    
    bool open(const char* host, int port, std::string& dbname, bool wrt);
    
    void close();

    int compareSchema(const std::string& schema);

    rdbModel::Rdb* m_rdb;
    bool m_write;
    bool m_dbg;
    /// optional database nickname
    std::string m_nick;  
  };
}

#endif
