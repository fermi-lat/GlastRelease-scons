/**
 * @class BuildEventStructure
 *
 * @brief Builds the McEventStructure TDS object from the existing base Monte Carlo 
 *        information output from G4Generator (McParticle and McPositionHit).
 *
 * @author The Tracking Software Group
 *
 * $Header$
 */
#include "GaudiKernel/SmartRefVector.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/IParticlePropertySvc.h"
#include "Event/MonteCarlo/McEventStructure.h"

#ifndef BuildEventStructure_h
#define BuildEventStructure_h

namespace Event {

class BuildEventStructure 
{
public:
    /// Standard Gaudi Tool interface constructor
    BuildEventStructure(IDataProviderSvc* dataSvc, IParticlePropertySvc* ppSvc);
   ~BuildEventStructure();

private:
    // Creates and fills a new McEventStructure object
    Event::McEventStructure* newMcEventStructure(IDataProviderSvc* dataSvc, IParticlePropertySvc* ppSvc);

    // This determines if input McParticle is a direct descendant of the primary
    //bool isPrimaryDaughter(const Event::McParticle* primary, const Event::McParticle* mcPart);

    // This finds all McParticles "associated" to the  input particle
    void findAssociated(Event::McEventStructure* mcEvent, const Event::McParticle* mcPart);

    // This "finds" the McParticle pointed at by mcPart in an McParticleRefVec
    bool find(const Event::McParticleRefVec::iterator begin, 
              const Event::McParticleRefVec::iterator end,
              const Event::McParticle*          mcPart);
};

};

#endif
