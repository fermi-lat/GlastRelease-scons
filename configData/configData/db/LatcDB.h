/** @class LatcDB
    @author M. Kocian
    

    Interface class to retrieve a filename from a database.

    $Id$
*/

#ifndef LATCDB_H
#define LATCDB_H

#include <string>

class LatcDB{
 public:
  /// Virtual destructor
  virtual ~LatcDB(){};
  /// virtual function to retrieve filename of type type given key key
  virtual const std::string getFilename(const char* type, unsigned key)=0;
};

#endif
