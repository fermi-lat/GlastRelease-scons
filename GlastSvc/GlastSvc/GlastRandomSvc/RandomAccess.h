#ifndef _GlastSvc_RandomAccess_H
#define _GlastSvc_RandomAccess_H 1

/** @class RandomAccess
 *
* @brief Gaudi Tool to set and get address of CLHEP random engines
* 
* This tool contains methods for 1) setting the CLHEP random engine 
* associated with GLAST shared libraries that use CLHEP random 
* number generators, and 2) returning the address of the random
* engine associated with a GLAST shared library. It is meant to
* be inherited by a tool (generally containing no additional methods)
* created in packages that generates shared libraries, one or more
* of whose methods use CLHEP random numbers. GlastRandomSvc, in
* it's initialize method, queries for the existence of tools that 
* inherit from RandomAccess and 1) sets the CLHEP engine type and 2) obtains 
* the address of the CLHEP random engine associated with the shared 
* library containing the tool. 
*     
*
* @authors Toby Burnett, Karl Young
*
* $Header
*/

// Include files
#include <string>

#include "GlastSvc/GlastRandomSvc/IRandomAccess.h"

#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/MsgStream.h"

#include "CLHEP/Random/Random.h"
#include "CLHEP/Random/DualRand.h"
#include "CLHEP/Random/TripleRand.h"
#include "CLHEP/Random/RandEngine.h"
#include "CLHEP/Random/DRand48Engine.h"
#include "CLHEP/Random/Hurd160Engine.h"
#include "CLHEP/Random/Hurd288Engine.h"
#include "CLHEP/Random/JamesRandom.h"
#include "CLHEP/Random/MTwistEngine.h"
#include "CLHEP/Random/RanecuEngine.h"
#include "CLHEP/Random/Ranlux64Engine.h"
#include "CLHEP/Random/RanluxEngine.h"
#include "CLHEP/Random/RanshiEngine.h"

class RandomAccess : public AlgTool, virtual public IRandomAccess 
{
 public:
  
  RandomAccess( const std::string& type, const std::string& name,
                const IInterface* parent) : AlgTool( type, name, parent ) 
    {
      // Declare additional interface
      declareInterface<IRandomAccess>(this);
    }
  
  ~RandomAccess()
    {
    }

  /// implement to define sources: will be called from FluxSvc, etc...
  HepRandomEngine* getTheEngine()
    {
      // Open the message log
      // MsgStream log( msgSvc(), name() );
      // log << MSG::DEBUG << "Getting address of random engine" << endreq;
      return HepRandom::getTheEngine();
    }

  void setTheEngine(std::string  engineName) 
    {  
      // Open the message log
      MsgStream log( msgSvc(), name() );
      int foundAnEngine = 0;
  
      if (engineName == std::string("DRand48Engine")) { 
	HepRandom::setTheEngine(new DRand48Engine);
	foundAnEngine = 1;
      }
      if (engineName == std::string("DualRand")) { 
	HepRandom::setTheEngine(new DualRand);
	foundAnEngine = 1;
      }
      if (engineName == std::string("HepJamesRandom")) { 
	HepRandom::setTheEngine(new HepJamesRandom);
	foundAnEngine = 1;
      }  
      if (engineName == std::string("Hurd160Engine")) { 
	HepRandom::setTheEngine(new Hurd160Engine);
	foundAnEngine = 1;
      }
      if (engineName == std::string("Hurd288Engine")) { 
	HepRandom::setTheEngine(new Hurd288Engine);
	foundAnEngine = 1;
      }
      if (engineName == std::string("MTwistEngine")) { 
	HepRandom::setTheEngine(new MTwistEngine);
	foundAnEngine = 1;
      }
      if (engineName == std::string("RandEngine")) { 
	HepRandom::setTheEngine(new RandEngine);
	foundAnEngine = 1;
      }
      if (engineName == std::string("RanecuEngine")) { 
	HepRandom::setTheEngine(new RanecuEngine);
	foundAnEngine = 1;
      }
      if (engineName == std::string("RanluxEngine")) { 
	HepRandom::setTheEngine(new RanluxEngine);
	foundAnEngine = 1;
      }
      if (engineName == std::string("Ranlux64Engine")) { 
	HepRandom::setTheEngine(new Ranlux64Engine);
	foundAnEngine = 1;
      }
      if (engineName == std::string("RanshiEngine")) { 
	HepRandom::setTheEngine(new RanshiEngine);
	foundAnEngine = 1;
      }
      if (engineName == std::string("TripleRand")) { 
	HepRandom::setTheEngine(new TripleRand);
	foundAnEngine = 1;
      }
      if (foundAnEngine == 0) { 
	log << MSG::ERROR << "Couldn't set random engine to "
	    << engineName 
	    << endreq;
      }
      else {
	// log << MSG::DEBUG << "Set random engine to "
	//    << engineName 
	//    << endreq;
      }
    }
};

#endif // _GlastSvc_RandomAccess_H
