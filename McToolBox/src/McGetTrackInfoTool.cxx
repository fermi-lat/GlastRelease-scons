/** 
 * @class McGetTrackInfoTool
 *
 * @brief A Gaudi tool for extracting information from the McParticle - TkrTrack and 
 *        TkrTrackHit relational tables for studying pattern reconstruction issues. 
 *
 *        NOTE: The interface for this Gaudi tool is defined 
 *              in GlastSvc/MonteCarlo/IMcGetTrackInfoTool
 *
 *        This version under construction!! 2/19/2004
 *
 * Created 1-Feb-2004
 * 
 * @author Tracy Usher
 *
 * $Header$
 */


#include "GlastSvc/MonteCarlo/IMcGetTrackInfoTool.h"
#include "GlastSvc/MonteCarlo/IMcGetEventInfoTool.h"

#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/GaudiException.h" 
#include "GaudiKernel/IParticlePropertySvc.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/ParticleProperty.h"
#include "GaudiKernel/AlgTool.h"

#include "Event/TopLevel/EventModel.h"
#include "Event/TopLevel/MCEvent.h"


class McGetTrackInfoTool : public AlgTool, virtual public IMcGetTrackInfoTool 
{
public:
    /// Standard Gaudi Tool interface constructor
    McGetTrackInfoTool(const std::string& type, const std::string& name, const IInterface* parent);
    virtual ~McGetTrackInfoTool() {}
	
    /// @brief Intialization of the tool
    StatusCode                    initialize();

    /// @brief Return the number of Monte Carlo tracks
    int                           getNMcParticles(const Event::TkrTrack*);

    /// @brief Return a vector of McParticles which have hits in the tracker
    const Event::McParticleRefVec getMcParticleVec(const Event::TkrTrack*);

    /// @brief Return TkrTrack <-> McParticle "best" match
    const Event::McParticle*      getBestMcParticle(const Event::TkrTrack*);

    /// @brief Return TkrTrack <-> McParticle "best" match
    const Event::TkrTrack*        getBestTkrTrack(const Event::McParticle*);

    /// @brief Return the number of hits on a given TkrTrack track associated to
    ///        a given McParticle
    int                           getNumMcHits(const Event::TkrTrack*, const Event::McParticle*);

private:
    /// Method for updating data
    const bool updateData();

    /// Pointer to the service which keeps track of the particle properties (most useful)
    IParticlePropertySvc*          m_ppsvc;

    /// Event Service member directly useable by concrete classes.
    IDataProviderSvc*              m_dataSvc;

    /// Use the McGetEventInfoTool for stuff that has already been calculated
    IMcGetEventInfoTool*           m_eventInfo;

    /// Event number to key on loading new tables
    TimeStamp                      m_time;           // Will use this when conversion to it complete
    int                            m_lastEventNo;    // backup for now

    /// Pointers to the Monte Carlo information for a single event
    Event::McPartToTkrTrackHitTab* m_partToTkrTrackHitTab;
    Event::McPartToTkrTrackTab*    m_partToTkrTrackTab;

    /// Null particle reference
    Event::McParticle              m_nullParticle;
};


//static ToolFactory<McGetTrackInfoTool> s_factory;
//const IToolFactory& McGetTrackInfoToolFactory = s_factory;
DECLARE_TOOL_FACTORY(McGetTrackInfoTool);

//
// Class constructor, no initialization here
//

McGetTrackInfoTool::McGetTrackInfoTool(const std::string& type, const std::string& name, const IInterface* parent) :
                       AlgTool(type, name, parent), m_time(0), m_lastEventNo(-1), m_nullParticle()
{
    //Declare additional interface
    declareInterface<IMcGetTrackInfoTool>(this);

    m_partToTkrTrackHitTab = 0;
    m_partToTkrTrackTab = 0;

	return;
}

//
// Initialization of the tool here
//

StatusCode McGetTrackInfoTool::initialize()
{	
    AlgTool::initialize();
    StatusCode sc   = StatusCode::SUCCESS;

    if( (sc = service("ParticlePropertySvc", m_ppsvc)).isFailure() ) 
    {
        throw GaudiException("Service [ParticlePropertySvc] not found", name(), sc);
    }

    if( (sc = service("EventDataSvc", m_dataSvc)).isFailure() ) 
    {
        throw GaudiException("Service [EventDataSvc] not found", name(), sc);
    }

    if ( (sc = toolSvc()->retrieveTool("McGetEventInfoTool", m_eventInfo)).isFailure() )
    {
        throw GaudiException("Tool [McGetEventInfoTool] not found", name(), sc);
    }

  return sc;
}

const bool McGetTrackInfoTool::updateData()
{
    // Assume success
    bool loaded = true;

    // Retrieve the pointer to the McEventStructure
    SmartDataPtr<Event::MCEvent> mcEvent(m_dataSvc,EventModel::MC::Event);

    if (mcEvent.ptr())
    {
        if (mcEvent->time() != m_time || mcEvent->getSequence() != m_lastEventNo)
        {
            m_time        = mcEvent->time();
            m_lastEventNo = mcEvent->getSequence();

            // Clean up the last table (if one)
            if (m_partToTkrTrackHitTab) delete m_partToTkrTrackHitTab;
            if (m_partToTkrTrackTab)    delete m_partToTkrTrackTab;

            // Retrieve the McParticle to Pat Cand Hit relational table
            SmartDataPtr<Event::McPartToTkrTrackHitTabList> TkrTrackHitTable(m_dataSvc,EventModel::MC::McPartToTkrTrackHitTab);
            m_partToTkrTrackHitTab = new Event::McPartToTkrTrackHitTab(TkrTrackHitTable);

            // Retrieve the McParticle to Pat Cand Track relational table
            SmartDataPtr<Event::McPartToTkrTrackTabList> trackTable(m_dataSvc,EventModel::MC::McPartToTkrTrackTab);
            m_partToTkrTrackTab = new Event::McPartToTkrTrackTab(trackTable);
        }
    }
    else
    {
        m_time = TimeStamp(0);
        loaded = false;
    }

    return loaded;
}

int McGetTrackInfoTool::getNMcParticles(const Event::TkrTrack* patCand)
{
    int numParticles = 0;

    if (updateData())
    {
        // Find the hits associated with this particle
        Event::McPartToTkrTrackVec partVec = m_partToTkrTrackTab->getRelBySecond(patCand);

        numParticles = partVec.size();
    }

    return numParticles;
}

const Event::McParticleRefVec McGetTrackInfoTool::getMcParticleVec(const Event::TkrTrack* patCand)
{
    Event::McParticleRefVec mcPartVec;

    mcPartVec.clear();

    if (updateData())
    {
        // Find the hits associated with this particle
        Event::McPartToTkrTrackVec partVec = m_partToTkrTrackTab->getRelBySecond(patCand);

        for(Event::McPartToTkrTrackVec::iterator partIter = partVec.begin(); partIter != partVec.end(); partIter++)
        {
            const Event::McParticle* mcPart = (*partIter)->getFirst();

            mcPartVec.push_back(mcPart);
        }
    }

    return mcPartVec;
}

const Event::McParticle* McGetTrackInfoTool::getBestMcParticle(const Event::TkrTrack* patCand)
{
    Event::McParticle* mcPartBest =  0;
    int                nBestHits  = -1;

    if (updateData())
    {
        // Find the hits associated with this particle
        Event::McPartToTkrTrackVec partVec = m_partToTkrTrackTab->getRelBySecond(patCand);

        for(Event::McPartToTkrTrackVec::iterator partIter = partVec.begin(); partIter != partVec.end(); partIter++)
        {
            Event::McPartToTkrTrackRel* rel         = *partIter;
            std::vector<std::string>      infoVec     = rel->getInfos();
            int                           infoVecSize = infoVec.size();

            if (infoVecSize > nBestHits)
            {
                mcPartBest = rel->getFirst();
                nBestHits  = infoVecSize;
            }
        }
    }

    return mcPartBest;
}

const Event::TkrTrack* McGetTrackInfoTool::getBestTkrTrack(const Event::McParticle* mcPart)
{
    const Event::TkrTrack* patCand   =  0;
    int                      nBestHits = -1;

    if (updateData())
    {
        // Find the hits associated with this particle
        Event::McPartToTkrTrackVec partVec = m_partToTkrTrackTab->getRelByFirst(mcPart);

        for(Event::McPartToTkrTrackVec::iterator partIter = partVec.begin(); partIter != partVec.end(); partIter++)
        {
            Event::McPartToTkrTrackRel* rel         = *partIter;
            std::vector<std::string>      infoVec     = rel->getInfos();
            int                           infoVecSize = infoVec.size();

            if (infoVecSize > nBestHits)
            {
                patCand    = rel->getSecond();
                nBestHits  = infoVecSize;
            }
        }
    }

    return patCand;
}

int McGetTrackInfoTool::getNumMcHits(const Event::TkrTrack* patCand, const Event::McParticle* mcPart)
{
    int numMcHits = 0;

    if (updateData())
    {
        // Find the hits associated with this particle
        Event::McPartToTkrTrackVec partVec = m_partToTkrTrackTab->getRelBySecond(patCand);

        for(Event::McPartToTkrTrackVec::iterator partIter = partVec.begin(); partIter != partVec.end(); partIter++)
        {
            Event::McPartToTkrTrackRel* rel         = *partIter;

            if (mcPart == rel->getFirst())
            {
                std::vector<std::string> infoVec = rel->getInfos();
                
                numMcHits = infoVec.size();

                break;
            }
        }
    }

    return numMcHits;
}
