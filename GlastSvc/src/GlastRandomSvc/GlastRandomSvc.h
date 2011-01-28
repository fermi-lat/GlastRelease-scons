/*
 * @class GlastRandomSvc
 *
 * @brief Gaudi Service for setting the random engines and seeds
 * for shared libraries that use random numbers (via CLHEP) 
 *
 * This service, in its initialize() method, collects the adresses 
 * of all tools that implement the  IRandomAccess interface (one in each
 * Dll that uses random numbers). The RandomAccess tool lives in 
 * GlastRandomSvc. The initialize() method also sets the random engine
 * for each of its Dll's either via the job options parameter
 * RandomEngine or the default which is currently TripleRand. The
 * handle() methods listens for BeginEvent events via the
 * IncidentListener service and increments the run and particle 
 * sequence numbers, sets those in the MCEvent header, then sets the
 * seed for each of the Dll's that use randoms, based on the run and
 * particle sequence numbers.
 * 
 *
 * @authors Toby Burnett, Karl Young, Michael Kuss
 *
 * $Header$
 */

#ifndef _GlastRandomSvc_H
#define _GlastRandomSvc_H 1

#include <vector>
#include <map>
#include "GaudiKernel/Service.h"
#include "GaudiKernel/IIncidentListener.h"
#include "GaudiKernel/Property.h"

#include "CLHEP/Random/Random.h"
#include "GlastRandomSeed.h"
#include "GlastSvc/GlastRandomSvc/IRandomAccess.h"

class IDataProviderSvc;;

class GlastRandomSvc : public Service,
    virtual public IIncidentListener
{
public:

    //! Create the service
    //! @param name The name of the service
    //!
    GlastRandomSvc ( const std::string& name, ISvcLocator* al );

    virtual ~GlastRandomSvc ();

    /// Handles incidents, implementing IIncidentListener interface
    void handle(const Incident& inc);

    /// required for Gaudi service
    StatusCode initialize ();

    /// required for Gaudi service
    StatusCode finalize ();

    /// action function that can be called by a (local) client-- will set all appropriate seeds
    /** @param run the current run number
    @param seq the sequence number
    */
    void applySeeds(int run, int seq); 

    /// singleton interface for local access
    static GlastRandomSvc * instance();

private:  
    static GlastRandomSvc* s_instance;

    CLHEP::HepRandomEngine* createEngine(std::string  engineName) ;

    /// Data members
    // store Engine names and addresses in a map
    typedef   std::map< std::string, CLHEP::HepRandomEngine* > EngineMap;
    EngineMap m_engineMap;
    StringProperty    m_randomEngine;
    std::string       m_runNumberString;
    IntegerProperty   m_runNumber;
    IntegerProperty   m_initialSequenceNumber;
    IntegerProperty   m_numberOfRuns;
    /// for setting the seed manually (run number seeded runs only)
    IntegerProperty   m_seed;
    int m_sequenceNumber;
    int m_maxEventId;
    unsigned int m_engineNum;

    // file name of seeds to be read in
    // if m_seedFile is empty, procceed to normal running of Gleam
    StringProperty m_seedFile;

    // Output seeds at the end of event, only used in testing whether event can
    // be regenerated
    StringProperty m_endSeedFile;
    std::ofstream m_output;

    // false: seed once at the beginning of a run
    // true: seed for each event
    BooleanProperty m_eventSeeded;

    // check if choice of m_runNumber, m_numberOfRuns, and
    // m_initialSequenceNumber ensure distinct random number sequences
    // false: issue warnings only
    // true: additionally, exit with StatusCode::FAILURE
    BooleanProperty m_pedantic;
    
    // seeds read from the file
    std::vector<GlastRandomSeed> m_seeds;

    IDataProviderSvc* m_eventSvc;

    //! list of pointers to setflag functions
    std::vector< IRandomAccess::SetFlag> m_setFlagPointers;

};

#endif // _GlastRandomSvc_H
