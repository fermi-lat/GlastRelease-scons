/** @file FluxObs.h
@brief definition of the class FluxObs

$Header$

*/
#ifndef _FluxObs_H
#define _FluxObs_H 1

#include "GaudiKernel/IToolSvc.h"
#include "GaudiKernel/Property.h"

#include "FluxSvc/IRegisterSource.h"
#include "FluxSvc/IFluxSvc.h"

//class FluxSvc;

/** @class FluxObs
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
class FluxObs : public IToolSvc::Observer
{
public:

    FluxObs();

    virtual ~FluxObs();

    void onCreate(IAlgTool& tool);
   
    void onRetrieve(IAlgTool& tool) { }

    void setFluxSvc(IFluxSvc *fluxSvc) { m_fluxSvc = fluxSvc; }

private:  

    IFluxSvc *m_fluxSvc;

};

#endif // _FluxObs_H

