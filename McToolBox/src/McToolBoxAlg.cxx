/** 
 * @class McToolBoxAlg
 *
 * @brief A Gaudi algorithm the drives the basic Monte Carlo event structure building and then 
 *        relates McParticles to their associated McPositionHits (for simple "MC" track finding
 *        at a later stage)
 *        If succesful, the TDS output of this algorithm will be an McEventStructure object (to
 *        characterize the event) and an McPartToPosHitTab relational table (relating the important
 *        McParticles to their associated McPositionHits). 
 *        **NOTE**: This algorithm will only operate in "pruneCal" or "full" McParticle pruning 
 *                  modes with fully reliable results only produced in "full" mode.
 *
 * Created 1-Feb-2004
 * 
 * @author Tracy Usher
 *
 * $Header$
 */

// Include standard gaudi stuff for an algorithm
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/DataObject.h"
#include "GaudiKernel/GaudiException.h" 

// Include stuff necessary to set up a display routine (if needed)
#include "gui/DisplayControl.h"
#include "gui/DisplayRep.h"
#include "GuiSvc/IGuiSvc.h"
#include "gui/GuiMgr.h"

// Include the tool interface for doing the real work...
#include "GlastSvc/MonteCarlo/IMcBuildRelTablesTool.h"
#include "GlastSvc/MonteCarlo/IMcGetEventInfoTool.h"

// Include our MC specific stuff
#include "Event/TopLevel/EventModel.h"
#include "Event/MonteCarlo/McEventStructure.h"
#include "Event/MonteCarlo/McRelTableDefs.h"
#include "Event/Recon/TkrRecon/TkrCluster.h"

//----------------------------------------------
//
//   MciTracksRep
//
//   Displays the Monte Carlo "tracks" (cluster to cluster)
//----------------------------------------------
//             Tracy Usher, SLAC, July 28, 2003
//----------------------------------------------
//##########################################################
class McTracksRep : public gui::DisplayRep
//##########################################################
{
public:
    //! Constructor of this form must be provided
    McTracksRep(IMcGetEventInfoTool* mcInfo);
    virtual ~McTracksRep() {}

    void update();

private:
    void drawTrack(const Event::McPartToClusPosHitVec hitVec, const std::string& color);

    IMcGetEventInfoTool* m_McInfo;
};

// Define the Algorithm
class McToolBoxAlg : public Algorithm
{
public:
    // Standard Gaudi Algorithm constructor format
    McToolBoxAlg(const std::string& name, ISvcLocator* pSvcLocator); 
    virtual ~McToolBoxAlg() {}

    // The thee phases in the life of a Gaudi Algorithm
    StatusCode initialize();
    StatusCode execute();
    StatusCode finalize();
    
private:
    // Pointer to the tool for building everything
    IMcBuildRelTablesTool* m_McInfo;
    IMcGetEventInfoTool*   m_McTracks;
};

//static const AlgFactory<McToolBoxAlg>  Factory;
//const IAlgFactory& McToolBoxAlgFactory = Factory;
DECLARE_ALGORITHM_FACTORY(McToolBoxAlg);

McToolBoxAlg::McToolBoxAlg(const std::string& name, ISvcLocator* pSvcLocator) :
Algorithm(name, pSvcLocator)  
{ 
}

StatusCode McToolBoxAlg::initialize()
{
    // Purpose and Method: Initialization method for the Tracker MC relations algorithm
    // Inputs:  None
    // Outputs:  StatusCode upon completetion
    // Dependencies: none
    // Restrictions and Caveats:  None

    StatusCode sc = StatusCode::SUCCESS;

    MsgStream  log(msgSvc(), name());

    setProperties();

    //Hunt for the tool for creating the relational tables
    sc = toolSvc()->retrieveTool("McBuildRelTablesTool", m_McInfo);

    //Hunt for the tool which returns mines the invaluable information from the tables
    sc = toolSvc()->retrieveTool("McGetEventInfoTool", m_McTracks);

    //Look for the gui service
    IGuiSvc* guiSvc = 0;

    sc = service("GuiSvc", guiSvc);
    if( sc.isSuccess() )  
    {
        gui::DisplayControl& display = guiSvc->guiMgr()->display();
        
        gui::DisplayControl::DisplaySubMenu& mcmenu = display.subMenu("McInfo");

        mcmenu.add(new McTracksRep(m_McTracks), "Monte Carlo Tracks");
    }
    else sc = StatusCode::SUCCESS;
   
    return sc;
}

StatusCode McToolBoxAlg::execute()
{
    // Purpose and Method: Called each event to build first the tables relating McParticles, McPositionHits
    //                     and Tkr Clusters to form MC tracks, then to relate this information to the 
    //                     Tracker Recon output.
    // Inputs:  None
    // Outputs:  Relations tables in the Tracker TDS, StatusCode upon completetion
    // Dependencies: None
    // Restrictions and Caveats:  None

    m_McInfo->buildEventStructure();
    m_McInfo->buildMonteCarloTracks();
//    m_McInfo->buildMcPatCandRelations();

    return StatusCode::SUCCESS;
}



StatusCode McToolBoxAlg::finalize()
{   
    return StatusCode::SUCCESS;
}

//------------------------------------------------------------------------------
/// Algorithm parameters which can be set at run time must be declared.
/// This should be done in the constructor.

//#############################################################################
McTracksRep::McTracksRep(IMcGetEventInfoTool* mcInfo)
//#############################################################################
{
    m_McInfo = mcInfo;
}
//-------------------- private ----------------------
//##############################################
void McTracksRep::update()
//##############################################
{
    // Recover pointer to the McEventStructure
    const Event::McEventStructure* mcEvent = m_McInfo->getMcEventStructure();

    // If mcEvent doesn't exist then there is no data to look at
    if (mcEvent)
    {

        // If the primary is charged we draw it first
        if (m_McInfo->getClassificationBits() & Event::McEventStructure::CHARGED) 
            drawTrack(m_McInfo->getMcPartTrack(mcEvent->getPrimaryParticle()),"purple");

        // Now draw the secondaries
        Event::McParticleRefVec::const_iterator partIter;

        for(partIter = mcEvent->beginSecondaries(); partIter != mcEvent->endSecondaries(); partIter++)
        {
            drawTrack(m_McInfo->getMcPartTrack(*partIter), "red");
        }

        for(partIter = mcEvent->beginAssociated(); partIter != mcEvent->endAssociated(); partIter++)
        {
            drawTrack(m_McInfo->getMcPartTrack(*partIter), "orange");
        }
    }

    return;
}

void McTracksRep::drawTrack(const Event::McPartToClusPosHitVec hitVec, const std::string& color)
{
    gui::DisplayRep* pDisplay = this;

    // If no hits then we aren't doing any drawing...
    if (hitVec.size() > 0)
    {
        // Recover the McParticle associated with this track
        const Event::McParticle* mcPart = (hitVec.front())->getFirst();

        // Set line segment start and end points
        Point x0(mcPart->initialPosition().x(),mcPart->initialPosition().y(),mcPart->initialPosition().z());

        Event::McPartToClusPosHitVec::const_iterator hitIter;
        for(hitIter = hitVec.begin(); hitIter != hitVec.end(); hitIter++)
        {
            Event::McPartToClusPosHitRel* mcHitRel = *hitIter;
            Event::McPositionHit*         mcPosHit = (mcHitRel->getSecond())->getSecond();
            const HepPoint3D              hitPos   = m_McInfo->getPosition(mcPosHit);

            Point  x1(hitPos.x(),hitPos.y(),hitPos.z());

            // do them in this order, so that the connection doesn't cover the track
        
            pDisplay->setColor(color);
            pDisplay->moveTo(x0);
            pDisplay->lineTo(x1); 

            x0 = x1;
        }
    }

    return;
}
