/** @class LatcDBImpl
    @author M. Kocian
    

    Implementation to retrieve a filename from MOOT 

    $Id$
*/
#ifndef LATCDBIMPL_H
#define LATCDBIMPL_H

#include "configData/db/LatcDB.h"
#include <string>
#include "mootCore/MootQuery.h"

using namespace MOOT;

class LatcDBImpl:public LatcDB{
 public:
  /// Default constructor
  LatcDBImpl();
  /// function to retrieve filename of type type given key key.
  const std::string getFilename(const char* type, unsigned key);
 private:
  MootQuery* m_mq;
};

#endif
