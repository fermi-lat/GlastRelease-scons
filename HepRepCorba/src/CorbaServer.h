// $Header:
//
// Description: 
//
// Author:R.Giannitrapni

#ifndef CORBASERVER_H
#define CORBASERVER_H

#include "HepRepSvc/IServer.h"

class IRegistry;
class ISvcAdapter;

/** 
 * @class IServer
 *
 * @brief A server interface for interactive HepRep implementations
 *    
 */
class CorbaServer: public IServer
{
 public:
  virtual void run();

 private:

};

#endif //CORBASERVER_H
