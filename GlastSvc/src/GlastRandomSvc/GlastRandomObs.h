/** @file GlastRandomObs.h
@brief definition of the class GlastRandomObs

$Header$

*/
#ifndef _GlastRandomObs_H
#define _GlastRandomObs_H 1

#include <vector>
#include <map>
#include "GaudiKernel/IToolSvc.h"
#include "GaudiKernel/Property.h"

#include "CLHEP/Random/Random.h"
#include "GlastRandomSeed.h"
#include "GlastSvc/GlastRandomSvc/IRandomAccess.h"


/** @class GlastRandomObs
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
* @authors Toby Burnett, Karl Young
*
* $Header$
*/
class GlastRandomObs : public IToolSvc::Observer
{
public:

    GlastRandomObs();

    virtual ~GlastRandomObs();

    void onCreate(IAlgTool& tool);
   
    void onRetrieve(IAlgTool& tool) { }

    const std::vector< IRandomAccess::SetFlag>& getSetFlagPtrs() const { return m_setFlagPointers; }
    
    typedef   std::map< std::string, CLHEP::HepRandomEngine* > EngineMap;

    const EngineMap& getEngineMap() const { return m_engineMap; }

    void setRandomEngine(const std::string engine) { m_randomEngine = engine; }

    void setRunNumber(int runNum) { m_runNumber = runNum; }

    void setSequenceNumber(int seqNum) { m_sequenceNumber = seqNum; }

    void setAutoSeed(bool autoSeed) { m_autoSeed = autoSeed; }

    int getEngineNum() const { return m_engineNum; }

    /// action function that can be called by a (local) client-- will set all appropriate seeds
    /** @param run the current run number
    @param seq the sequence number
    */
    //void applySeeds(int run, int seq); 


private:  

    CLHEP::HepRandomEngine* createEngine(std::string  engineName) ;

    /// Data members
    EngineMap m_engineMap;
    std::string    m_randomEngine;

    int m_runNumber;
    int m_sequenceNumber;
    int m_engineNum;

    bool m_autoSeed;
 
    //! list of pointers to setflag functions
    std::vector< IRandomAccess::SetFlag> m_setFlagPointers;

};

#endif // _GlastRandomObs_H

