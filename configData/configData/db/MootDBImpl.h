/** @class MootDBImpl
    @author M. Kocian
    

    Implementation to retrieve sets of files from MOOT

    $Id$
*/
#ifndef MOOTDBIMPL_H
#define MOOTDBIMPL_H

#include "configData/db/LatcDB.h"
#include <string>
#include "mootCore/MootQuery.h"
#include "mootCore/FileDescrip.h"

using namespace MOOT;

class MootDBImpl:public LatcDB{
 public:
  /// Default constructor
  MootDBImpl();

  /// function to retrieve filename of type type given key key.
  const std::string getFilename(const char* type, unsigned key);

  /// function to retrieve the list of filenames given a MOOT key
  const std::vector<std::string> getFilenameList(unsigned key);

 private:
  MootQuery* m_mq;
};

#endif
