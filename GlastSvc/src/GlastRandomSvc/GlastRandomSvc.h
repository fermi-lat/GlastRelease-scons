
#ifndef _GlastRandomSvc_H
#define _GlastRandomSvc_H 1

#include <map>
#include "GaudiKernel/Service.h"
#include "GaudiKernel/IEvtSelector.h"
#include "GaudiKernel/IIncidentListener.h"
#include "CLHEP/Random/Random.h"

#include "GlastSvc/GlastRandomSvc/IGlastRandomSvc.h"

/** @class GlastRandomSvc
* @brief This is a Gaudi Service for setting the random engines
* and seeds for shared libraries that use random numbers (via CLHEP) 
* 
* @author Toby Burnett, Karl Young
* $Header
*/



class GlastRandomSvc : public Service,
virtual public IIncidentListener,
virtual public IGlastRandomSvc
{
public:
  
  /// Handles incidents, implementing IIncidentListener interface
  void handle(const Incident& inc);
  
  //! Create the service
  //! @param name The name of the service
  //!
  GlastRandomSvc ( const std::string& name, ISvcLocator* al );
  
  virtual ~GlastRandomSvc ();
  
  /// queryInterface - for implementing a Service this is necessary
  StatusCode queryInterface(const IID& riid, void** ppvUnknown);
  
  /// return the service type
  const IID& type() const;
  
  StatusCode initialize ();
  
  StatusCode finalize ();
  
  // set random seed for engine
  StatusCode setSeeds();

  // query for available tools and make the map
  StatusCode getEngines();

 private:  
  /// Data members
  // store Engine names and addresses in a map
  std::map< std::string, HepRandomEngine* > m_engineMap;
  std::string     m_randomEngine;
};

#endif // _GlastRandomSvc_H


