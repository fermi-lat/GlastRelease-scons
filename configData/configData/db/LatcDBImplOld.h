/** @class LatcDBImplOld
    @author M. Kocian
    

    Implementation to retrieve a filename from MOOT (using an old function that will be replaced soon.

    $Id$
*/
#ifndef LATCDBIMPLOLD_H
#define LATCDBIMPLOLD_H

#include "configData/db/LatcDB.h"
#include <string>
#include "mootCore/MootQuery.h"

using namespace MOOT;

class LatcDBImplOld:public LatcDB{
 public:
  /// Default constructor
  LatcDBImplOld();
  /// function to retrieve filename of type type given key key.
  const std::string getFilename(const char* type, unsigned key);
  /// function to retrieve the list of filenames
  const std::vector<std::string> getFilenameList(unsigned key);
 private:
  MootQuery* m_mq;
};

#endif
