// $Header$
#ifndef MOOT_MootSys_h
#define MOOT_MootSys_h

#include <string>

/** @file MootSys.h
    @author J. Bogart
  
    MOOT needs to do various things via system calls
    Might also handle other parts of building a config; TBD.
    This is a private class.  Do not expose via Python interface.
*/

namespace MOOT {
  class MootSys {
  public:
    /**
       Compute crc and length.
       returns:  true if successful
       This version expects a file path.  Would be preferable to
       take arbitrary stream as input.
    */
    static bool computeCRC(const std::string& fname, std::string &crc, 
                           std::string &sz);

    /** 
        Issue a command to the system
           @a cmd       Command to issue (input)
           @a out       Results of command which would normally go to std::out
           @a status    If command executed, its return.  Else -1.
           @a clear     If true, clear output string at start
    */
    static bool sysCmd(const std::string& cmd, std::string& out,
                       int& status, bool clear=true);

    /// Check whether supplied path corresponds to directory.
    /// Normally @a path will be absolute
    static bool isDir(const std::string& path);

    /// Make child directory of parent.  Normally @a path will be absolute
    static bool makeDir(const std::string& parent, const std::string& dirname);

    /// Expand environment variables (of form $(env) ) and append just
    /// portion of name to the right of last "/" to @a base.
    static void basename(const std::string& fullpath, std::string& base);

    /// Extract value of parameter from a string.  
    /// @arg   str   string to search
    /// @arg   pname   name of parameter whose value is sought
    /// @arg   pvalue  string passed by ref. to hold return value
    /// @arg   conn    string connecting parm name to value, defaults to "="
    /// @arg   sep     string separating one parm from next, defaults to " "  
    /// @return true if parm found
    static bool getParmValue(const std::string& str, const std::string& pname, 
                             std::string& pvalue, 
                             const std::string& conn=std::string("="),
                             const std::string& sep=std::string(" "));
  };

}
#endif
