#ifndef _GlastRandomSvc_H
#define _GlastRandomSvc_H 1

/** @class GlastRandomSvc
 *
* @brief Gaudi Service for setting the random engines and seeds
* for shared libraries that use random numbers (via CLHEP) 
* 
* This service, in its initialize() method, collects the adresses 
* of all tools that inherit form the RandomAccess tool (one in each
* Dll that uses random numbers). The RandomAccess tool lives in 
* GlastRandomAvc. The initialize() method also sets the random engine
* for each of its Dll's either via the job options parameter
* RandomEngine or the default which is currently TripleRand. The
* handle() methods listens for BeginEvent events via the
* IncidentListener service and increments the run and particle 
* sequence numbers, sets those in the MCEvent header, then sets the
* seed for each of the Dll's that use randoms, based on the run and
* particle sequence numbers.
* 
*
* @authors Toby Burnett, Karl Young
*
* $Header
*/

#include <iostream>
#include <vector>
#include <map>
#include "GaudiKernel/Service.h"
#include "GaudiKernel/IEvtSelector.h"
#include "GaudiKernel/IIncidentListener.h"
#include "GaudiKernel/Property.h"

#include "CLHEP/Random/Random.h"

#include "GlastSvc/GlastRandomSvc/IGlastRandomSvc.h"

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
  
  /// required for Gaudi service
  StatusCode initialize ();
  
  /// required for Gaudi service
  StatusCode finalize ();
  
 private:  
  /// Data members
  // store Engine names and addresses in a map
  std::map< std::string, HepRandomEngine* > m_engineMap;
  StringProperty    m_randomEngine;
  IntegerProperty   m_RunNumber;
  IntegerProperty   m_InitialSequenceNumber;
  int m_SequenceNumber;

  // file name of seeds to be read in
  // if m_seedFile is empty, procceed to normal running of Gleam
  string m_seedFile;

  struct Seed {
    long m_run;
    long m_seqNo;
  };
  friend istream& operator >> (istream& in, Seed& seed);
  friend ostream& operator << (ostream& in, const Seed& seed);

  // seeds read from the file
  std::vector<Seed> m_seeds;
};

extern istream& operator >> (istream& in, GlastRandomSvc::Seed& seed);
extern ostream& operator << (ostream& in, const GlastRandomSvc::Seed& seed);

#endif // _GlastRandomSvc_H


