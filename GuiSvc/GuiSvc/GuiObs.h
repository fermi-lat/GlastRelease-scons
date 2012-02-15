/* Definition of the class GuiObs

$Header$

*/
#ifndef _GuiObs_H
#define _GuiObs_H 1

#include "GaudiKernel/IToolSvc.h"
#include "GaudiKernel/Property.h"

#include "GuiSvc/IGuiTool.h"

/** @class GuiObs
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
class GuiObs : public IToolSvc::Observer
{
public:

    GuiObs();

    virtual ~GuiObs();

    void onCreate(const IAlgTool* tool);
   
    void onRetrieve(const IAlgTool* tool) { }

    void setGuiMgr(gui::GuiMgr* guiMgr) { m_guiMgr = guiMgr; }

private:  
    gui::GuiMgr* m_guiMgr;

};

#endif // _GuiObs_H

