// $Header:
//
// Description: 
//
// Author:R.Giannitrapni

#ifndef XMLSTREAMER_H
#define XMLSTREAMER_H

#include <string>
#include "HepRepSvc/IStreamer.h"

class XmlStreamer: public IStreamer
{
 public:
  /// The pure virtual method to save the HepRep to a filename
  virtual void saveHepRep(std::string);

};

#endif //XMLSTREAMER_H
