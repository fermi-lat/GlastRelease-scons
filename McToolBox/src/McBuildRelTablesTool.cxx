/** 
 * @class IMcBuildRelTablesTool
 *
 * @brief A Gaudi tool that drives the basic Monte Carlo event structure building and then 
 *        relates McParticles to their associated McPositionHits (for simple "MC" track finding
 *        at a later stage)
 *        If succesful, the TDS output of this tool will be an McEventStructure object (to
 *        characterize the event) and an McPartToPosHitTab relational table (relating the important
 *        McParticles to their associated McPositionHits). 
 *        **NOTE**: This tool will only operate in "pruneCal" or "full" McParticle pruning 
 *                  modes with fully reliable results only produced in "full" mode.
 *
 * Created 1-Feb-2004
 * 
 * @author Tracy Usher
 *
 * $Header$
 */

// Include standard gaudi stuff for an algorithm
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/GaudiException.h" 
#include "GaudiKernel/IParticlePropertySvc.h"
#include "GaudiKernel/ParticleProperty.h"
#include "GaudiKernel/AlgTool.h"

// Include our interface
#include "GlastSvc/MonteCarlo/IMcBuildRelTablesTool.h"
#include "Event/TopLevel/EventModel.h"

// Include the individual class definitions for doing the work
#include "BuildEventStructure.h"
#include "BuildMcTracks.h"
#include "BuildPatCandTab.h"

// Define the Tool
class McBuildRelTablesTool : public AlgTool, virtual public IMcBuildRelTablesTool 
{
public:
    // Standard Gaudi Algorithm constructor format
    McBuildRelTablesTool(const std::string& type, const std::string& name, const IInterface* parent); 
    virtual ~McBuildRelTablesTool() {}
	
    /// @brief Intialization of the tool
    StatusCode initialize();

    /// @Builds the event structure
    void buildEventStructure();

    /// @Builds tables relating Tracker "hits" to McParticles
    void buildMonteCarloTracks();

    /// @Builds tables relating MC tracks to TKR Pattern Candidate tracks
    void buildMcPatCandRelations();
    
private:

    /// Pointer to the service which keeps track of the particle properties (most useful)
    IParticlePropertySvc*    m_ppsvc;

    /// Event Service member directly useable by concrete classes.
    IDataProviderSvc*        m_dataSvc;
};

//static ToolFactory<McBuildRelTablesTool> s_factory;
//const IToolFactory& McBuildRelTablesToolFactory = s_factory;
DECLARE_TOOL_FACTORY(McBuildRelTablesTool);

McBuildRelTablesTool::McBuildRelTablesTool(const std::string& type, const std::string& name, const IInterface* parent) :
                       AlgTool(type, name, parent)
{ 
    //Declare additional interface
    declareInterface<IMcBuildRelTablesTool>(this);
}

StatusCode McBuildRelTablesTool::initialize()
{
    // Purpose and Method: Initialization method for the Tracker MC relations algorithm
    // Inputs:  None
    // Outputs:  StatusCode upon completetion
    // Dependencies: none
    // Restrictions and Caveats:  None

    StatusCode sc = StatusCode::SUCCESS;

    MsgStream log(msgSvc(), name());

    setProperties();

    if( (sc = service("ParticlePropertySvc", m_ppsvc)).isFailure() ) {
        throw GaudiException("Service [ParticlePropertySvc] not found", name(), sc);
    }

    if( (sc = service("EventDataSvc", m_dataSvc)).isFailure() ) {
        throw GaudiException("Service [EventDataSvc] not found", name(), sc);
    }
   
    return sc;
}

void McBuildRelTablesTool::buildEventStructure()
{
    // Purpose and Method: Called each event to build first the tables relating McParticles, McPositionHits
    //                     and Mci Clusters to form MC tracks, then to relate this information to the 
    //                     Tracker Recon output.
    // Inputs:  None
    // Outputs:  Relations tables in the Tracker TDS, StatusCode upon completetion
    // Dependencies: None
    // Restrictions and Caveats:  None

    // Always believe in success
    StatusCode sc = StatusCode::SUCCESS;

    Event::BuildEventStructure eventStructure(m_dataSvc, m_ppsvc);

    // finished
    return;
}

void McBuildRelTablesTool::buildMonteCarloTracks()
{
    // Purpose and Method: Called each event to build first the tables relating McParticles, McPositionHits
    //                     and Mci Clusters to form MC tracks, then to relate this information to the 
    //                     Tracker Recon output.
    // Inputs:  None
    // Outputs:  Relations tables in the Tracker TDS, StatusCode upon completetion
    // Dependencies: None
    // Restrictions and Caveats:  None

    // Always believe in success
    StatusCode sc = StatusCode::SUCCESS;

    Event::BuildMcTracks mcTracks(m_dataSvc);

    // finished
    return;
}

void McBuildRelTablesTool::buildMcPatCandRelations()
{
    // Purpose and Method: Called each event to build first the tables relating McParticles, McPositionHits
    //                     and Mci Clusters to form MC tracks, then to relate this information to the 
    //                     Tracker Recon output.
    // Inputs:  None
    // Outputs:  Relations tables in the Tracker TDS, StatusCode upon completetion
    // Dependencies: None
    // Restrictions and Caveats:  None

    // Always believe in success
    StatusCode sc = StatusCode::SUCCESS;

    Event::BuildPatCandTab mcPatCandRel(m_dataSvc);

    // finished
    return;
}

