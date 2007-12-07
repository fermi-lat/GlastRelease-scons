/** @class LatcDB
    @author M. Kocian
    

    Interface class to retrieve a filename from a database.

    $Id$
*/

#ifndef LATCDB_H
#define LATCDB_H

#include <string>
#include <vector>

class LatcDB{
 public:
  /// Virtual destructor
  virtual ~LatcDB(){};
  /// virtual function to retrieve filename of type type given key key
  virtual const std::string getFilename(const char* type, unsigned key)=0;
  /// function to retrieve the list of filenames
  virtual const std::vector<std::string> getFilenameList(unsigned key)=0;
};

#endif
