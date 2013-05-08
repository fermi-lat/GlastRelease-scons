// $Header$

// Try out utilities in MootSys
// from .my.cnf

#include <iostream>
#include <string>
#include "mootCore/MootSys.h"
#include "facilities/commonUtilities.h"


int main(int, char**)    {

  MOOT::MootSys  myUtil;
  std::string out;
  std::string cmd("pwd");
  int status;

  facilities::commonUtilities::setupEnvironment();

  // Issue a command that returns a single line
  bool ret = myUtil.sysCmd(cmd, out, status);

  std::cout << "Command '" << cmd << "'  completed with status " << status
            << " and returned this output: " << std::endl;
  std::cout << out;
  std::cout << "sysCmd returned " << ret << "." << std::endl << std::endl;
  std::cout.flush();

  cmd = std::string("printenv PATH");
  ret = myUtil.sysCmd(cmd, out, status);

  std::cout << "Command '" << cmd << "'  completed with status " << status
            << " and returned this output: " << std::endl;
  std::cout << out;
  std::cout << "sysCmd returned " << ret << "." << std::endl << std::endl;
  std::cout.flush();

  cmd = std::string("which zcompress");
  ret = myUtil.sysCmd(cmd, out, status);
  std::cout << "Command '" << cmd << "'  completed with status " << status
            << " and returned this output: " << std::endl;
  std::cout << out;
  std::cout << "sysCmd returned " << ret << "." << std::endl << std::endl;
  std::cout.flush();


  // Try again with multi-line output
  /*   this will have to go in MOOT, not mootCore, for meaningful output
  cmd = std::string("fmx directory /mordor");
  ret = myUtil.sysCmd(cmd, out, status);
  std::cout << "Command '" << cmd << "'  completed with status " << status
            << " and returned this output: " << std::endl;
  std::cout << out;
  std::cout << "sysCmd returned " << ret << "." << std::endl << std::endl;
  std::cout.flush();

  //  cmd = std::string("/u/ey/jrb/dev500/mootdev/FMX/V2-0-1/cmt/fmx directory /mordor");

  cmd = std::string("fmx badcmd /mordor");
  ret = myUtil.sysCmd(cmd, out, status);
  std::cout << "Command '" << cmd << "'  completed with status " << status
            << " and returned this output: " << std::endl;
  std::cout << out;
  std::cout << "sysCmd returned " << ret << "." << std::endl << std::endl;
  std::cout.flush();


  cmd = std::string("fmx directory /nosuchinstrument");
  ret = myUtil.sysCmd(cmd, out, status);
  std::cout << "Command '" << cmd << "'  completed with status " << status
            << " and returned this output: " << std::endl;
  std::cout << out;
  std::cout.flush();

  */
  std::string crc, sz;

  ret = myUtil.computeCRC("../cmt/require17ments", crc, sz);

  if (ret) {
    std::cout << "Command 'cksum requirements' returned crc=" << crc
              << " and size =" << sz << std::endl;
  }
  else std::cout << "computeCRC failed " << std::endl;
  std::cout.flush();

  // Try to make a subdirectory
  bool ok = myUtil.makeDir("/u/ey/jrb/dev500", "junkdir");   // should succeed
  std::cout << "Attempt to make /u/ey/jrb/dev500/junkdir returned status "
            << ok << std::endl;

  ok = myUtil.makeDir ("/u/ey/jrb", "transfer");  // already exists
  std::cout << "Attempt to make /u/ey/jrb/transfer returned status "
            << ok << std::endl;

  ok = myUtil.makeDir("/u/ey/jrb/xfr", "new");   // should fail
  std::cout << "Attempt to make /u/ey/jrb/xfr/new returned status "
            << ok << std::endl;

  std::string dbtype;
  std::string master;
  std::string local;
  std::string cmdOut;
  /*  Move to mootCore for meaningful response
  int cmdStatus;
  ok = MOOT::MootSys::sysCmd("fmx show connection", cmdOut, cmdStatus);
  if (ok) {
    ok = MOOT::MootSys::getParmValue(cmdOut, std::string("Database"), dbtype);
    if (ok) std::cout << "fmx database is: " << dbtype << std::endl;
    else std::cout << "fmx database not found" <<  std::endl;
    ok = MOOT::MootSys::getParmValue(cmdOut, std::string("master"), master);
    if (ok) std::cout << "fmx master is: " << master << std::endl;
    else std::cout << "fmx master not found" <<  std::endl;
    ok = MOOT::MootSys::getParmValue(cmdOut, std::string("local"), local);
    if (ok) std::cout << "fmx local node is: " << local << std::endl;
    else std::cout << "fmx local node not found" <<  std::endl;
  }
  */
  cmdOut = "...with junk=boo key=74)";
  std::string keyStr, junk;

  ok = MOOT::MootSys::getParmValue(cmdOut, std::string("junk"), junk,
                                   "=", " ),");
  if (ok) std::cout << "junk is: " << junk << std::endl;
  else std::cout << "junk not found" <<  std::endl;

  ok = MOOT::MootSys::getParmValue(cmdOut, std::string("key"), keyStr,
                                   "=", " ),");
  if (ok) std::cout << "key is: " << keyStr << std::endl;
  else std::cout << "key not found" <<  std::endl;

  return 0;
}
