/** @class LatcDBImplFile
    @author M. Kocian
    

    Implementation to retrieve a filename that was set previously to fit in with the DB interface

    $Id$
*/
#ifndef LATCDBIMPLFILE_H
#define LATCDBIMPLFILE_H

#include "configData/db/LatcDB.h"
#include <string>


class LatcDBImplFile:public LatcDB{
 public:
  /// Default constructor
  LatcDBImplFile();
  /// function to retrieve filename of type type given key key.
  const std::string getFilename(const char* type, unsigned key);
  /// function to set the filename
  void setFilename(const std::string filename);
 private:
  std::string m_filename;
};

#endif
