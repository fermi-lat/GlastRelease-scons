// $Header$
#include <cstdio>
//#include <cstdlib>
#include <iostream>
#include <vector>
#include "mootCore/MootSys.h"
#include "facilities/Util.h"
#include <sys/stat.h>
#include <cerrno>

namespace MOOT {

  bool MootSys::computeCRC(const std::string& fname, 
                           std::string &crc, std::string &sz) {
#ifndef WIN32
    std::string cmd("cksum ");
    cmd = cmd + fname;

    std::string out;
    int status;

    if (!sysCmd(cmd, out, status)) return false;

    if (status != 0) { // something went wrong
      std::cerr << "Command " << cmd << " returned with status " 
                << status << std::endl;
      return false;
    }

    // Otherwise command executed successfully; need to parse output
    std::vector<std::string> tokens;
    
    facilities::Util  myFacUtil;
    
    myFacUtil.stringTokenize(out, " ", tokens);

    // should return crc, size, and filename

    if (tokens.size() != 3) {
      std::string msg2("computeCRC: bad output for file ");
      msg2 = msg2 + fname + "\n";
      std::cerr << msg2;
      return 0;
    }

    crc = tokens[0];
    sz = tokens[1];
    return true;
#else
    return false;           // shouldn't call from Windows
#endif
  }

  bool MootSys::sysCmd(const std::string& cmd, std::string& out, 
                       int& status, bool clear) {
#ifndef WIN32
    FILE *stream;

    status = -1;   // initialized for failure
    if (clear) out.clear();
    stream = popen(cmd.c_str(), "r");

    if (!stream) {
      std::string msg("Failed popen for command: ");
      msg = msg + cmd + "\n";
      std::cerr << msg;
      return false;
    }
    const int BUFSIZE = 256;
    char buffer[BUFSIZE];

    char* rd;
    bool notDone=true;
    while (notDone) {
      rd = fgets(buffer, BUFSIZE, stream);
      if (!rd) {
        notDone = false;
        if (ferror(stream)) {
          std::string msg("Failed to read output for command: ");
          msg = msg + cmd + "\n";
          std::cerr << msg;
          pclose(stream);
          return false;
        }
      } else  { // got some output
        out += std::string(buffer);
      }
    }
    status = pclose(stream);
    return (status == 0);
#else
    return false;         // shouldn't call from Windows
#endif
  }

  // Return true if file exists and is a directory; else false
  bool MootSys::isDir(const std::string& path) {
#ifndef WIN32
    struct stat buf;

    int ret = stat(path.c_str(), &buf);
    if (ret != 0) return false;           // couldn't stat the file

    return S_ISDIR(buf.st_mode);
#else
    return false;     // shouldn't call from Windows. Put out error??
#endif
  }

  bool MootSys::makeDir(const std::string& parent, 
                        const std::string& dirname) {

    // As long as the archive is in afs space only user mode bits
    // are used.
    //
    // Can set up acls to allow Insert, but not Delete or Write
    // maybe change owner so that original owner can't delete it.

    // First see if it already exists.  If so, just return true.
#ifndef WIN32
    std::string path = parent + std::string("/") + dirname;
    if (isDir(path)) return true;

    if (!(isDir(parent))) return false;

    // need only user permissions in afs space but for devel. on nfs
    // give everyone read access
    static mode_t mode = S_IRWXU | S_IRWXG | S_IRWXO;   

    int ret = mkdir(path.c_str(), mode);
    int errCode = errno;
    if (ret) {  
      std::cerr << "MootSys::makeDir  System call mkdir failed with code "
                << errCode << std::endl;
    }
    return (ret == 0);
#else
    return false;      // shouldn't call from Windows
#endif
  }

  void MootSys::basename(const std::string& fullPath, std::string& base)
  {
#ifndef WIN32
    unsigned slashPos = fullPath.rfind('/');
#else
    unsigned slashPos = fullPath.rfind('\\');
#endif
    base += fullPath.substr(slashPos + 1);
    return;
  }

  bool MootSys::getParmValue(const std::string& str, const std::string& pname, 
                             std::string& pvalue, const std::string& conn,
                             const std::string& sep) {
    std::string search = pname + conn;
    pvalue.clear();
    //      std::string search("Database=");
    unsigned start = str.find(search);
    if (start >= str.size()) return false;

    start += search.size();
    unsigned end = str.find_first_of(sep, start);
    if (end >= str.size()) end = str.size();
    pvalue.assign(str, start, end - start);
    return true;
  }

}

