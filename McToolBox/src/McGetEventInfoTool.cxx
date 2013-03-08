/** 
 * @class McGetEventInfoTool
 *
 * @brief A Gaudi tool for extracting information from the McParticle - McPositionHit -
 *        TkrCluster - etc. relational tables. The primary aim is to aid in studies of 
 *        track reconstruction performance, particular pattern recognition.  
 *
 *        NOTE: The interface for this Gaudi tool is defined 
 *              in GlastSvc/MonteCarlo/IMcGetEventInfoTool
 *
 *        This version under construction!! 2/19/2004
 *
 * Created 1-Feb-2004
 * 
 * @author Tracy Usher
 *
 * $Header$
 */


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
#include "Event/MonteCarlo/McEventStructure.h"
#include "Event/Recon/TkrRecon/TkrCluster.h"


class McGetEventInfoTool : public AlgTool, virtual public IMcGetEventInfoTool 
{
public:
    /// Standard Gaudi Tool interface constructor
    McGetEventInfoTool(const std::string& type, const std::string& name, const IInterface* parent);
    virtual ~McGetEventInfoTool() {}
	
    /// @brief Intialization of the tool
    StatusCode                  initialize();

    /// @brief Following methods return information on type of event and particles
    /// @brief Also available directly from the TDS
    const Event::McEventStructure* getMcEventStructure();
    /// @brief Returns the number of Monte Carlo tracks in the tracker
    int                         getNumMcTracks();
    /// @brief Return a vector of McParticles which have hits in the tracker
    const Event::McParticleRefVec getMcTrackVector();
    /// @brief Returns information about the event
    const unsigned long         getClassificationBits();
    /// @brief Returns primary McParticle
    const Event::McParticleRef  getPrimaryParticle();
    /// @brief Returns secondary particles
    int                         getNumSecondaries();
    const Event::McParticleRef  getSecondary(int mcPartIdx);
    /// @brief Returns associated particles
    int                         getNumAssociated();
    const Event::McParticleRef  getAssociated(int mcPartIdx);

    /// @brief Returns a vector of hits associated as one McParticle track
    const Event::McPartToClusPosHitVec getMcPartTrack(const Event::McParticleRef mcPart);

    /// @brief Returns the layer number of a given McPositionHit on the track
    const int                   getTrackHitLayer(const Event::McParticleRef mcPart, int hitIdx);

    /// @brief Following methods return information about specific tracks
    /// @brief Returns number of Tracker (cluster) hits for a given track
    const int                   getNumClusterHits(const Event::McParticleRef mcPart);
    /// @brief Returns number of shared Tracker (cluster) hits for a given track
    const int                   getNumSharedTrackHits(const Event::McParticleRef mcPart);
    /// @brief Compares two tracks and returns information on shared hits (if any)
    const int                   getNumGaps(const Event::McParticleRef mcPart);
    const int                   getGapSize(const Event::McParticleRef mcPart, int gapIdx);
    const int                   getGapStartHitNo(const Event::McParticleRef mcPart, int gapIdx);
    /// @brief Returns the "straightness" of a given track
    const double                getTrackStraightness(const Event::McParticleRef mcPart, int firstHitIdx=0, int lastHitIdx=40);
    /// @brief Returns the "direction" defined by the first set of hits on a track
    const CLHEP::Hep3Vector     getTrackDirection(const Event::McParticleRef mcPart, int firstHitIdx=0, int lastHitIdx=40);
    /// @brief Returns track energy loss information (in the tracker only)
    const double                getTrackTotEneLoss(const Event::McParticleRef mcPart);
    const double                getTrackELastHit(const Event::McParticleRef mcPart);
    const double                getTrackBremELoss(const Event::McParticleRef mcPart, int& nTotRad);
    const double                getTrackDeltaELoss(const Event::McParticleRef mcPart, int& nTotDlta, int& nHitDlta);
    const int                   getTrackDeltaRange(const Event::McParticleRef mcPart, double& aveRange, double& maxRange);
    /// @brief Compares two tracks and returns information on shared hits (if any)
    const unsigned int          getSharedHitInfo(const Event::McParticleRef mcPart);
    const unsigned int          getSharedHitInfo(const Event::McParticleRef mcPart1, const Event::McParticleRef mcPart2);

    /// @brief Calculate position info from McPositionHits related to a given cluster
    const HepPoint3D getPosition(const Event::TkrCluster* cluster);
    const HepPoint3D getPosition(const Event::McPositionHit* mcPosHit);

private:
    /// Method for updating data
    const bool updateData();

    /// Pointer to the service which keeps track of the particle properties (most useful)
    IParticlePropertySvc*         m_ppsvc;

    /// Event Service member directly useable by concrete classes.
    IDataProviderSvc*             m_dataSvc;

    /// Event number to key on loading new tables
    TimeStamp                     m_time;           // Will use this when conversion to it complete
    int                           m_lastEventNo;    // backup for now

    /// Pointers to the Monte Carlo information for a single event
    Event::McEventStructure*      m_mcEvent;
    Event::ClusMcPosHitTab*       m_clusHitTab;
    Event::McPartToClusTab*       m_partClusTab;
    Event::McPartToClusPosHitTab* m_partHitTab;

    /// Null particle reference
    Event::McParticle             m_nullParticle;
};


//static ToolFactory<McGetEventInfoTool> s_factory;
//const IToolFactory& McGetEventInfoToolFactory = s_factory;
DECLARE_TOOL_FACTORY(McGetEventInfoTool);

//
// Class constructor, no initialization here
//

McGetEventInfoTool::McGetEventInfoTool(const std::string& type, const std::string& name, const IInterface* parent) :
                       AlgTool(type, name, parent), m_time(0), m_lastEventNo(-1), m_nullParticle()
{
    //Declare additional interface
    declareInterface<IMcGetEventInfoTool>(this);

    m_mcEvent     = 0;
    m_clusHitTab  = 0;
    m_partClusTab = 0;
    m_partHitTab  = 0;

	return;
}

//
// Initialization of the tool here
//

StatusCode McGetEventInfoTool::initialize()
{	
  AlgTool::initialize();
  StatusCode sc   = StatusCode::SUCCESS;

  if( (sc = service("ParticlePropertySvc", m_ppsvc)).isFailure() ) {
        throw GaudiException("Service [ParticlePropertySvc] not found", name(), sc);
  }

  if( (sc = service("EventDataSvc", m_dataSvc)).isFailure() ) {
        throw GaudiException("Service [EventDataSvc] not found", name(), sc);
  }

  return sc;
}

const bool McGetEventInfoTool::updateData()
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

            // Retrieve the pointer to the McEventStructure
            m_mcEvent = SmartDataPtr<Event::McEventStructure>(m_dataSvc,EventModel::MC::McEventStructure);

            // Clean up the last table (if one)
            if (m_partHitTab)  delete m_partHitTab;
            if (m_partClusTab) delete m_partClusTab;
            if (m_clusHitTab)  delete m_clusHitTab;

            // Retrieve the TkrCluster <-> McPositionHit table
            SmartDataPtr<Event::ClusMcPosHitTabList> clusTable(m_dataSvc,EventModel::Digi::TkrClusterHitTab);
            m_clusHitTab = new Event::ClusMcPosHitTab(clusTable);

            // Retrieve the McParticle <-> TkrCluster table
            SmartDataPtr<Event::McPartToClusTabList> hitTable(m_dataSvc,EventModel::MC::McPartToClusTab);
            m_partClusTab = new Event::McPartToClusTab(hitTable);

            // Retrieve the McParticle to hit relational table
            SmartDataPtr<Event::McPartToClusPosHitTabList> partTable(m_dataSvc,EventModel::MC::McPartToClusHitTab);
            m_partHitTab = new Event::McPartToClusPosHitTab(partTable);
        }
    }
    else
    {
        m_time = TimeStamp(0);
        loaded = false;
    }

    return loaded;
}


//
// Define a small class which can be used by the std::sort algorithm 
//
class CompareTrackHits 
{
public:
    bool operator()(Event::McPartToClusPosHitRel *left, Event::McPartToClusPosHitRel *right)
    {
        bool leftTest   = false;

        // Extract the TkrCluster <-> McPositionHit relation 
        const Event::ClusMcPosHitRel* mcHitLeft  = left->getSecond();
        const Event::ClusMcPosHitRel* mcHitRight = right->getSecond();

        // Extract the McPositionHit embedded in this relation
        const Event::McPositionHit*   mcPosHitLeft  = mcHitLeft->getSecond();
        const Event::McPositionHit*   mcPosHitRight = mcHitRight->getSecond();

        // If McPositionHits found, sort is by the particle's time of flight
        if (mcPosHitLeft && mcPosHitRight)
        {
            leftTest = mcPosHitLeft->timeOfFlight() < mcPosHitRight->timeOfFlight();
        }

        return leftTest;
    }
private:
};

//
// Return pointer to the McEventStructure
//
const Event::McEventStructure* McGetEventInfoTool::getMcEventStructure()
{
    Event::McEventStructure* mcEvent = 0;

    if (updateData()) mcEvent = m_mcEvent;

    return mcEvent;
}

const Event::McParticleRefVec McGetEventInfoTool::getMcTrackVector()
{
    Event::McParticleRefVec hitVec;
    hitVec.clear();

    if (updateData())
    {
        hitVec = m_mcEvent->getTrackVector();
    }

    return hitVec;
}

//
// How many Monte Carlo tracks in the event?
//

int McGetEventInfoTool::getNumMcTracks()
{
    // Always believe in success
    StatusCode sc = StatusCode::SUCCESS;

    // By default, no tracks
    int numMcTracks = 0;

    // If it doesn't exist then we need to build the MC structure
    if (updateData())
    {
        // Obtain the vector of McParticle references from McEventStructure
        Event::McParticleRefVec                 trackVec = m_mcEvent->getTrackVector();
        Event::McParticleRefVec::const_iterator trackVecIter;

        for(trackVecIter = trackVec.begin(); trackVecIter != trackVec.end(); trackVecIter++)
        {
            // Find the hits associated with this particle
            Event::McPartToClusPosHitVec hitVec = m_partHitTab->getRelByFirst(*trackVecIter);

            // Don't bother if really too few hits
            if (hitVec.size() > 4) numMcTracks++;
        }

        /*
        // If primary particle is charged then count if it is a track
        if (m_mcEvent->getClassificationBits() & Event::McEventStructure::CHARGED) 
        {
            // Find the hits associated with this particle
            Event::McPartToClusPosHitVec hitVec = m_partHitTab->getRelByFirst(m_mcEvent->getPrimaryParticle());

            // Don't bother if really too few hits
            if (hitVec.size() > 4) numMcTracks++;
        }

        // Now look at the secondaries
        Event::McParticleRefVec::const_iterator partIter;

        for(partIter = m_mcEvent->beginSecondaries(); partIter != m_mcEvent->endSecondaries(); partIter++)
        {
            // Find the hits associated with this particle
            Event::McPartToClusPosHitVec hitVec = m_partHitTab->getRelByFirst(*partIter);

            // Don't bother if really too few hits
            if (hitVec.size() > 4) numMcTracks++;
        }

        // Finally, any associated tracks
        for(partIter = m_mcEvent->beginAssociated(); partIter != m_mcEvent->endAssociated(); partIter++)
        {
            // Find the hits associated with this particle
            Event::McPartToClusPosHitVec hitVec = m_partHitTab->getRelByFirst(*partIter);

            // Don't bother if really too few hits
            if (hitVec.size() > 4) numMcTracks++;
        }
        */
    }

    return numMcTracks;
}

const unsigned long McGetEventInfoTool::getClassificationBits()
{
    if (updateData())
    {
        return m_mcEvent->getClassificationBits();
    }

    return 0;
}

const Event::McParticleRef  McGetEventInfoTool::getPrimaryParticle()
{
    if (updateData())
    {
        return m_mcEvent->getPrimaryParticle();
    }

    return &m_nullParticle;
}

int McGetEventInfoTool::getNumSecondaries()
{
    if (updateData()) return m_mcEvent->getNumSecondaries();

    return 0;
}

const Event::McParticleRef  McGetEventInfoTool::getSecondary(int mcPartIdx)
{
    if (updateData())
    {
        Event::McParticleRefVec::const_iterator refVec = m_mcEvent->beginSecondaries();

        if (mcPartIdx >= 0 && mcPartIdx < m_mcEvent->getNumSecondaries()) return refVec[mcPartIdx];
    }

    return &m_nullParticle;
}

int McGetEventInfoTool::getNumAssociated()
{
    if (updateData()) return m_mcEvent->getNumAssociated();

    return 0;
}

const Event::McParticleRef  McGetEventInfoTool::getAssociated(int mcPartIdx)
{
    if (updateData())
    {
        Event::McParticleRefVec::const_iterator refVec = m_mcEvent->beginAssociated();

        if (mcPartIdx >= 0 && mcPartIdx < m_mcEvent->getNumAssociated()) return refVec[mcPartIdx];
    }

    return &m_nullParticle;
}


const Event::McPartToClusPosHitVec McGetEventInfoTool::getMcPartTrack(const Event::McParticleRef mcPart)
{
    Event::McPartToClusPosHitVec hitVec;
    hitVec.clear();

    if (updateData())
    {
        hitVec = m_partHitTab->getRelByFirst(mcPart);
        std::sort(hitVec.begin(),hitVec.end(),CompareTrackHits());
    }

    return hitVec;
}

const int McGetEventInfoTool::getNumClusterHits(const Event::McParticleRef mcPart)
{
    int numTrackHits = 0;

    if (updateData())
    {
        Event::McPartToClusPosHitVec trackVec = m_partHitTab->getRelByFirst(mcPart);
        std::sort(trackVec.begin(),trackVec.end(),CompareTrackHits());

        Event::McPartToClusPosHitVec::const_iterator trackVecIter;

        for(trackVecIter = trackVec.begin(); trackVecIter != trackVec.end(); trackVecIter++)
        {
            Event::TkrCluster* cluster = ((*trackVecIter)->getSecond())->getFirst();

            if (cluster) numTrackHits++;
        }
    }

    return numTrackHits;
}

const int McGetEventInfoTool::getNumSharedTrackHits(const Event::McParticleRef mcPart)
{
    int numSharedHits = 0;

    if (updateData())
    {
        Event::McPartToClusPosHitVec trackVec = m_partHitTab->getRelByFirst(mcPart);
        std::sort(trackVec.begin(),trackVec.end(),CompareTrackHits());

        Event::McPartToClusPosHitVec::const_iterator trackVecIter;

        for(trackVecIter = trackVec.begin(); trackVecIter != trackVec.end(); trackVecIter++)
        {
            Event::TkrCluster*     cluster = ((*trackVecIter)->getSecond())->getFirst();
            Event::McPartToClusVec partVec = m_partClusTab->getRelBySecond(cluster);

            if (partVec.size() > 1) numSharedHits++;
        }
    }

    return numSharedHits;
}

const int McGetEventInfoTool::getNumGaps(const Event::McParticleRef mcPart)
{
    int numGaps = 0;

    if (updateData())
    {
        Event::McPartToClusPosHitVec trackVec = m_partHitTab->getRelByFirst(mcPart);

        // No need to proceed if not enough hits track
        if (trackVec.size() > 1)
        {
            std::sort(trackVec.begin(),trackVec.end(),CompareTrackHits());

            Event::McPartToClusPosHitVec::const_iterator trackVecIter = trackVec.begin();
            Event::McPositionHit*                        layerHit     = ((*trackVecIter++)->getSecond())->getSecond();
            idents::VolumeIdentifier                     volId        = layerHit->volumeID();
            int                                          lastLayer    = 2*volId[4] - 1 + volId[6];

            for(; trackVecIter != trackVec.end(); trackVecIter++)
            {
                layerHit = ((*trackVecIter)->getSecond())->getSecond();
                volId    = layerHit->volumeID();

                int tkrLayer = 2*volId[4] - 1 + volId[6];

                if (abs(tkrLayer - lastLayer) > 1) numGaps++;

                lastLayer = tkrLayer;
            }
        }
    }

    return numGaps;
}

const int McGetEventInfoTool::getGapSize(const Event::McParticleRef mcPart, int gapIdx)
{
    int gapSize = 0;

    if (updateData())
    {
        Event::McPartToClusPosHitVec trackVec = m_partHitTab->getRelByFirst(mcPart);

        // No need to proceed if not enough hits track
        if (trackVec.size() > 1)
        {
            std::sort(trackVec.begin(),trackVec.end(),CompareTrackHits());

            Event::McPartToClusPosHitVec::const_iterator trackVecIter = trackVec.begin();
            Event::McPositionHit*                        layerHit     = ((*trackVecIter++)->getSecond())->getSecond();
            idents::VolumeIdentifier                     volId        = layerHit->volumeID();
            int                                          lastLayer    = 2*volId[4] - 1 + volId[6];
            int                                          gapNum       = 0;

            for(; trackVecIter != trackVec.end(); trackVecIter++)
            {
                layerHit = ((*trackVecIter)->getSecond())->getSecond();
                volId    = layerHit->volumeID();

                int tkrLayer = 2*volId[4] - 1 + volId[6];

                if (abs(tkrLayer - lastLayer) > 1)
                {
                    if (gapNum++ == gapIdx) gapSize = abs(tkrLayer - lastLayer) - 1;
                }

                lastLayer = tkrLayer;
            }
        }
    }

    return gapSize;
}

const int McGetEventInfoTool::getGapStartHitNo(const Event::McParticleRef mcPart, int gapIdx)
{
    int gapStartHitNo = 0;

    if (updateData())
    {
        Event::McPartToClusPosHitVec trackVec = m_partHitTab->getRelByFirst(mcPart);

        // No need to proceed if not enough hits track
        if (trackVec.size() > 1)
        {
            std::sort(trackVec.begin(),trackVec.end(),CompareTrackHits());

            Event::McPartToClusPosHitVec::const_iterator trackVecIter = trackVec.begin();
            Event::McPositionHit*                        layerHit     = ((*trackVecIter++)->getSecond())->getSecond();
            idents::VolumeIdentifier                     volId        = layerHit->volumeID();
            int                                          lastLayer    = 2*volId[4] - 1 + volId[6];
            int                                          gapNum       = 0;
            int                                          hitNo        = 0;

            for(; trackVecIter != trackVec.end(); trackVecIter++)
            {
                layerHit = ((*trackVecIter)->getSecond())->getSecond();
                volId    = layerHit->volumeID();

                int tkrLayer = 2*volId[4] - 1 + volId[6];

                if (abs(tkrLayer - lastLayer) > 1)
                {
                    if (gapNum++ == gapIdx) gapStartHitNo = hitNo + 1;
                }

                lastLayer = tkrLayer;
                hitNo++;
            }
        }
    }

    return gapStartHitNo;
}


const double McGetEventInfoTool::getTrackStraightness(const Event::McParticleRef mcPart, int firstHitIdx, int lastHitIdx)
{
    double trackAngle  = 0.0;
    int    numClusHits = 0;

    if (updateData())
    {
        Event::McPartToClusPosHitVec trackVec = m_partHitTab->getRelByFirst(mcPart);
        int                          numHits  = trackVec.size();

        // Only keep going if we have enough hits to calculate an angle
        if (numHits > 3 && lastHitIdx - firstHitIdx > 3)
        {
            // Sort the hits to go from start to end of track
            std::sort(trackVec.begin(),trackVec.end(),CompareTrackHits());

            // Iterators over the hits in the track
            Event::McPartToClusPosHitVec::const_iterator trackVecIter = firstHitIdx < numHits
                                                                      ? trackVec.begin() + firstHitIdx
                                                                      : trackVec.begin();
            Event::McPartToClusPosHitVec::const_iterator trackVecStop = lastHitIdx < numHits
                                                                      ? trackVec.begin() + lastHitIdx 
                                                                      : trackVec.end();

            // Retrieve the first TkrCluster hit and get its position
            Event::TkrCluster* cluster  = ((*trackVecIter++)->getSecond())->getFirst();

            // It is possible for an McPositionHit to not make a TkrCluster, watch for this
            while(!cluster && trackVecIter != trackVecStop)
            {
                cluster = ((*trackVecIter++)->getSecond())->getFirst();
            }

            // Retrieve the next TkrCluster hit, but be sure we have no over run the iterator
            Event::TkrCluster* nextClus = 0;
            while(!cluster && trackVecIter != trackVecStop)
            {
                nextClus  = ((*trackVecIter++)->getSecond())->getFirst();
            }

            // If valid cluster hits then proceed
            if (cluster && nextClus)
            {
                // Determine the hit position - based on the cluster but from the McPositionHits
                HepPoint3D hitLast = getPosition(cluster);
                HepPoint3D hitPos  = getPosition(nextClus);

                // Form a vector between these hits
                CLHEP::Hep3Vector vecLast = CLHEP::Hep3Vector(hitPos - hitLast).unit();

                // Set up for looping over the remaining hits
                hitLast = hitPos;

                numClusHits = 2;

                // Now loop over the rest of the hits
                for( ; trackVecIter != trackVecStop; trackVecIter++)
                {
                    // Get next cluster hit
                    cluster = ((*trackVecIter)->getSecond())->getFirst();

                    if (!cluster) continue;

                    numClusHits++;

                    // Retrieve the McSiLayerHit position
                    hitPos  = getPosition(cluster);

                    // New vector to current hit
                    CLHEP::Hep3Vector vecNext = CLHEP::Hep3Vector(hitPos - hitLast).unit();

                    // update the track angle
                    double vecAngle = vecNext.angle(vecLast);
                    trackAngle += vecAngle * vecAngle;

                    // update for next loop
                    hitLast = hitPos;
                    vecLast = vecNext;
                }

                // Ok, if we found more than 2 clusters then calculate the rms
                if (numClusHits > 2)
                {
                    double divisor = numClusHits - 2;

                    trackAngle /= divisor;
                    trackAngle  = sqrt(trackAngle);
                }
                else trackAngle = 0.;
            }
        }
    }

    return trackAngle;
}

const HepPoint3D McGetEventInfoTool::getPosition(const Event::McPositionHit* mcPosHit)
{
    double x = 0.;
    double y = 0.;
    double z = 0.;

    if (updateData())
    {
        x = 0.5* (mcPosHit->globalEntryPoint().x() + mcPosHit->globalExitPoint().x());
        y = 0.5* (mcPosHit->globalEntryPoint().y() + mcPosHit->globalExitPoint().y());
        z = 0.5* (mcPosHit->globalEntryPoint().z() + mcPosHit->globalExitPoint().z());
    }

    return HepPoint3D(x,y,z);
}


const HepPoint3D McGetEventInfoTool::getPosition(const Event::TkrCluster* cluster)
{
    double x = 0.;
    double y = 0.;
    double z = 0.;

    if (updateData())
    {
        Event::ClusMcPosHitVec hitVec = m_clusHitTab->getRelByFirst(cluster);
        Event::ClusMcPosHitVec::const_iterator hitVecIter = hitVec.begin();

        // If only one McPositionHit associated with the cluster then task is easy
        if (hitVec.size() == 1)
        {
            Event::McPositionHit* posHit = (*hitVecIter)->getSecond();

            // Cluster is from a noise hit can't happen here? Check anyway
            if (posHit)
            {
                x = 0.5* (posHit->globalEntryPoint().x() + posHit->globalExitPoint().x());
                y = 0.5* (posHit->globalEntryPoint().y() + posHit->globalExitPoint().y());
                z = 0.5* (posHit->globalEntryPoint().z() + posHit->globalExitPoint().z());
            }
        }
        else if (hitVec.size() > 1)
        {
            double xLow  =  100000.;
            double xHigh = -100000.;
            double yLow  =  100000.;
            double yHigh = -100000.;
            double zLow  =  100000.;
            double zHigh = -100000.;

            for( ;hitVecIter != hitVec.end(); hitVecIter++)
            {
                Event::McPositionHit* posHit = (*hitVecIter)->getSecond();

                // It is possible that a cluster could contain a noise hit (hence, the
                // McPositionHit in the relation could be null
                if (!posHit) continue;

                HepPoint3D entryPoint = posHit->globalEntryPoint();
                if (entryPoint.x() < xLow)  xLow  = entryPoint.x();
                if (entryPoint.x() > xHigh) xHigh = entryPoint.x();
                if (entryPoint.y() < yLow)  yLow  = entryPoint.y();
                if (entryPoint.y() > yHigh) yHigh = entryPoint.y();
                if (entryPoint.z() < zLow)  zLow  = entryPoint.z();
                if (entryPoint.z() > zHigh) zHigh = entryPoint.z();

                HepPoint3D exitPoint  = posHit->globalExitPoint();
                if (exitPoint.x() < xLow)   xLow  = exitPoint.x();
                if (exitPoint.x() > xHigh)  xHigh = exitPoint.x();
                if (exitPoint.y() < yLow)   yLow  = exitPoint.y();
                if (exitPoint.y() > yHigh)  yHigh = exitPoint.y();
                if (exitPoint.z() < zLow)   zLow  = exitPoint.z();
                if (exitPoint.z() > zHigh)  zHigh = exitPoint.z();
            }

            x = 0.5 * (xLow + xHigh);
            y = 0.5 * (yLow + yHigh);
            z = 0.5 * (zLow + zHigh);
        }
    }

    return HepPoint3D(x,y,z);
}


const CLHEP::Hep3Vector McGetEventInfoTool::getTrackDirection(const Event::McParticleRef mcPart, int firstHitIdx, int lastHitIdx)
{
    double xSlope = 0.;
    double ySlope = 0.;

    if (updateData())
    {
        Event::McPartToClusPosHitVec trackVec = m_partHitTab->getRelByFirst(mcPart);
        int                   numHits  = trackVec.size();

        // Only keep going if we have enough hits to calculate an angle
        if (numHits > 3 && lastHitIdx - firstHitIdx > 3)
        {
            std::sort(trackVec.begin(),trackVec.end(),CompareTrackHits());

            // Iterator over the hits in the track
            Event::McPartToClusPosHitVec::const_iterator trackVecIter = firstHitIdx < numHits
                                                               ? trackVec.begin() + firstHitIdx
                                                               : trackVec.begin();
            Event::McPartToClusPosHitVec::const_iterator trackVecStop = lastHitIdx < numHits
                                                               ? trackVec.begin() + lastHitIdx 
                                                               : trackVec.end();

            double xVals[2];
            double yVals[2];
            double zVals_x[2];
            double zVals_y[2];
            int    nHitsX = 0;
            int    nHitsY = 0;

            // Now loop over the rest of the hits
            for( ; trackVecIter != trackVecStop; trackVecIter++)
            {
                // Get the cluster
                const Event::TkrCluster* tkrClus  = ((*trackVecIter)->getSecond())->getFirst();

                // Watch out for no cluster!
                if (!tkrClus) continue;

                // Check to see which orientation we have, store info accordingly
                if (tkrClus->getTkrId().getView() == idents::TkrId::eMeasureX && nHitsX < 2)
                {
                    xVals[nHitsX]   = tkrClus->position().x();
                    zVals_x[nHitsX] = tkrClus->position().z();
                    if (nHitsX == 1 && zVals_x[0] == zVals_x[1]) break;
                    nHitsX++;
                }
                else if (nHitsY < 2)
                {
                    yVals[nHitsY]   = tkrClus->position().y();
                    zVals_y[nHitsY] = tkrClus->position().z();
                    if (nHitsY == 1 && zVals_y[0] == zVals_y[1]) break;
                    nHitsY++;
                }

                // No use looping forever if done
                if (nHitsX > 1 && nHitsY > 1) break;
            }

            // if enough info then calculate new slopes
            if (nHitsX == 2 && nHitsY == 2)
            {
                xSlope = (xVals[1] - xVals[0]) / (zVals_x[1] - zVals_x[0]);
                ySlope = (yVals[1] - yVals[0]) / (zVals_y[1] - zVals_y[0]);
            }
        }
    }

    return CLHEP::Hep3Vector(-xSlope,-ySlope,-1.).unit();
}

const int McGetEventInfoTool::getTrackHitLayer(const Event::McParticleRef mcPart, int hitIdx )
{
    int layer = 50;

    if (updateData())
    {
        Event::McPartToClusPosHitVec trackVec = m_partHitTab->getRelByFirst(mcPart);
        int                   numHits  = trackVec.size();

        // Only keep going if we have enough hits to calculate an angle
        if (numHits > 0)
        {
            std::sort(trackVec.begin(),trackVec.end(),CompareTrackHits());

            if (hitIdx < 0)        hitIdx = 0;
            if (hitIdx >= numHits) hitIdx = numHits - 1;

            Event::McPositionHit*  layerHit = (trackVec[hitIdx]->getSecond())->getSecond();

            const idents::VolumeIdentifier volId = layerHit->volumeID();

            if (volId[0] == 0 && volId[3] == 1 && volId.size() > 6)
            {
                int trayNum = volId[4];
                int botTop  = volId[6];

                layer  = 2 * trayNum + botTop - 1;
            }
        }
    }

    return layer;
}


const double McGetEventInfoTool::getTrackTotEneLoss(const Event::McParticleRef mcPart)
{
    // Return the total energy loss from particle creation to last hit in tracker
    double totEloss = mcPart->initialFourMomentum().e() - getTrackELastHit(mcPart);

    return totEloss;
}

const double McGetEventInfoTool::getTrackELastHit(const Event::McParticleRef mcPart)
{
    double                       ELastHit = mcPart->initialFourMomentum().e();
    Event::McPartToClusPosHitVec hitVec   = getMcPartTrack(mcPart);
    int                          hitSize  = hitVec.size();

    // If we have McPositionHits associated with this particle, return the energy at
    // the exit point of the last McPositionHit
    if (hitSize > 0)
    {
        // Sort the vector into the proper track order (time ordered)
        std::sort(hitVec.begin(),hitVec.end(),CompareTrackHits());

        // Now get the energy at the last tracker hit
        Event::McPartToClusPosHitRel* relLast  = hitVec.back();
        Event::McPositionHit*         lastHit  = (relLast->getSecond())->getSecond();

        ELastHit = lastHit->particleEnergy();


        //following unecessary?
        /*
        const SmartRefVector<Event::McPositionHit>* mcPosHitVec = lyrLast->getMcPositionHitsVec();

        // This loop to take into account that McPositionHits can be broken across a layer
        // Need to find the last one (lowest energy)
        SmartRefVector<Event::McPositionHit>::const_iterator posHitIter;
        for(posHitIter = mcPosHitVec->begin(); posHitIter != mcPosHitVec->end(); posHitIter++)
        {
            const Event::McPositionHit* posHit = *posHitIter;

            if (posHit->particleEnergy() < ELastHit) ELastHit = posHit->particleEnergy();
        }
        */
    }
    else
    {
        // We should not be getting to this point
        int j=0;
    }

    return ELastHit;
}

const double McGetEventInfoTool::getTrackBremELoss(const Event::McParticleRef mcPart, int& nTotRad)
{
    double radEloss = 0.;

    nTotRad = 0;

    // McPositionHits for looking at where the last hit is...
    Event::McPartToClusPosHitVec hitVec  = getMcPartTrack(mcPart);
    const Event::McPositionHit*  posHit  = 0;
    int                          hitSize = hitVec.size();

    // If there exist McPositionHits for this mcPart then proceed...
    if (hitSize > 0)
    {
        // Sort the vector into the proper track order
        std::sort(hitVec.begin(),hitVec.end(),CompareTrackHits());

        // Get the last Layer Hit from which we can extract the last McPositionHit
        Event::McPartToClusPosHitRel* relLast = hitVec.back();
        Event::McPositionHit*         posHit  = (relLast->getSecond())->getSecond();

        // Get the daughter vector
        const SmartRefVector<Event::McParticle>& daughterVec = mcPart->daughterList();

        // Loop over particle's daughters
        SmartRefVector<Event::McParticle>::const_iterator daughterVecIter;
        for(daughterVecIter = daughterVec.begin(); daughterVecIter != daughterVec.end(); daughterVecIter++)
        { 
            const Event::McParticle* daughter = *daughterVecIter;

            // Retrieve the volume identifier so we can check daughter particle origin volume
            idents::VolumeIdentifier iniVolId = const_cast<Event::McParticle*>(daughter)->getInitialId();

            // If the particle starts in the tracker then record the energy it takes away
            if (iniVolId[0] == 0 && iniVolId[3] == 1 && iniVolId.size() > 3 && 
                daughter->initialPosition().z() > posHit->globalExitPoint().z())
            {
                // Looking only gammas
                ParticleProperty* ppty = m_ppsvc->findByStdHepID( daughter->particleProperty() );

                // why do we have to check this? 
                if (!ppty) continue;

                if (ppty->particle() == "gamma")
                {
                    radEloss += daughter->initialFourMomentum().e();
                    nTotRad++;
                }
            }
        }
    }

    return radEloss;
}

const double McGetEventInfoTool::getTrackDeltaELoss(const Event::McParticleRef mcPart, int& nTotDlta, int& nHitDlta)
{
    double radEloss = 0.;

    nTotDlta = 0;
    nHitDlta = 0;

    // McPositionHits for looking at where the last hit is...
    Event::McPartToClusPosHitVec       hitVec  = getMcPartTrack(mcPart);
    const Event::McPositionHit* posHit  = 0;
    int                         hitSize = hitVec.size();

    // If there exist McPositionHits for this mcPart then proceed...
    if (hitSize > 0)
    {
        // Sort the vector into the proper track order
        std::sort(hitVec.begin(),hitVec.end(),CompareTrackHits());

        // Get the last Layer Hit from which we can extract the last McPositionHit
        Event::McPartToClusPosHitRel* relLast = hitVec.back();
        Event::McPositionHit*         posHit  = (relLast->getSecond())->getSecond();

        const SmartRefVector<Event::McParticle>& daughterVec = mcPart->daughterList();
        SmartRefVector<Event::McParticle>::const_iterator daughterVecIter;

        // Loop over particle's daughters
        for(daughterVecIter = daughterVec.begin(); daughterVecIter != daughterVec.end(); daughterVecIter++)
        { 
            const Event::McParticle* daughter = *daughterVecIter;

            idents::VolumeIdentifier iniVolId = const_cast<Event::McParticle*>(daughter)->getInitialId();

            // If the particle starts in the tracker then record the energy it takes away
            if (iniVolId[0] == 0 && iniVolId[3] == 1 && iniVolId.size() > 3 && 
                daughter->initialPosition().z() > posHit->globalExitPoint().z())
            {
                // Get the particle property
                ParticleProperty* ppty = m_ppsvc->findByStdHepID( daughter->particleProperty() );

                // why do we have to check this? 
                if (!ppty) continue;

                // Looking for energy carried away by delta rays...
                if (ppty->particle() != "gamma")
                {
                    radEloss += (daughter->initialFourMomentum().e() - ppty->mass());

                    nTotDlta++;
                    if (daughter->statusFlags() & Event::McParticle::POSHIT) nHitDlta++;
                }
            }
        }
    }

    return radEloss;
}

const int McGetEventInfoTool::getTrackDeltaRange(const Event::McParticleRef mcPart, double& aveRange, double& maxRange)
{
    int nTotHits = 0;

    aveRange = 0.;
    maxRange = 0.;

    // McPositionHits for looking at where the last hit is...
    Event::McPartToClusPosHitVec       hitVec  = getMcPartTrack(mcPart);
    const Event::McPositionHit* posHit  = 0;
    int                         hitSize = hitVec.size();

    // If there exist McPositionHits for this mcPart then proceed...
    if (hitSize > 0)
    {
        // Sort the vector into the proper track order
        std::sort(hitVec.begin(),hitVec.end(),CompareTrackHits());

        // Get the last Layer Hit from which we can extract the last McPositionHit
        Event::McPartToClusPosHitRel* relLast = hitVec.back();
        Event::McPositionHit*         posHit  = (relLast->getSecond())->getSecond();

        const SmartRefVector<Event::McParticle>& daughterVec = mcPart->daughterList();
        SmartRefVector<Event::McParticle>::const_iterator daughterVecIter;

        // Loop over particle's daughters
        for(daughterVecIter = daughterVec.begin(); daughterVecIter != daughterVec.end(); daughterVecIter++)
        { 
            const Event::McParticle* daughter = *daughterVecIter;

            idents::VolumeIdentifier iniVolId = const_cast<Event::McParticle*>(daughter)->getInitialId();

            // If the particle starts in the tracker then record the energy it takes away
            if (iniVolId[0] == 0 && iniVolId[3] == 1 && iniVolId.size() > 3 && 
                daughter->initialPosition().z() > posHit->globalExitPoint().z())
            {
                // Get the particle property
                ParticleProperty* ppty = m_ppsvc->findByStdHepID( daughter->particleProperty() );

                // why do we have to check this? 
                if (!ppty) continue;

                // Looking for energy carried away by delta rays...
                if (ppty->particle() != "gamma")
                {
                    HepPoint3D        iniPosition = daughter->initialPosition();
                    HepPoint3D        finPosition = daughter->finalPosition();
                    CLHEP::Hep3Vector traject     = finPosition - iniPosition;
                    double            range       = traject.mag();

                    aveRange += range;

                    if (range > maxRange) maxRange = range;

                    nTotHits++;
                }
            }
        }

        // Calculate the average range
        if (nTotHits > 0) aveRange = aveRange / nTotHits;
    }

    return nTotHits;
}

const unsigned int McGetEventInfoTool::getSharedHitInfo(const Event::McParticleRef mcPart)
{
    unsigned int hitInfo = 0;

    if (updateData())
    {
        Event::McPartToClusPosHitVec trackVec = m_partHitTab->getRelByFirst(mcPart);
        std::sort(trackVec.begin(),trackVec.end(),CompareTrackHits());

        Event::McPartToClusPosHitVec::const_iterator trackVecIter;
        int                                   bitIndex     = 0;
        int                                   numShared    = 0;

        for(trackVecIter = trackVec.begin(); trackVecIter != trackVec.end(); trackVecIter++)
        {
            Event::TkrCluster*     cluster = ((*trackVecIter)->getSecond())->getFirst();
            Event::McPartToClusVec partVec = m_partClusTab->getRelBySecond(cluster);

            if (partVec.size() > 1) 
            {
                if (numShared < 15) numShared++;
                else                hitInfo |= 0x08000000;
                if (bitIndex < 24)  hitInfo |= 1 << bitIndex;
                else                hitInfo |= 0x00800000;
            }

            bitIndex++;
        }

        hitInfo = (hitInfo << 4) + numShared;
    }

    return hitInfo;
}


const unsigned int McGetEventInfoTool::getSharedHitInfo(const Event::McParticleRef mcPart1, const Event::McParticleRef mcPart2)
{
    unsigned int hitInfo = 0;

    if (updateData())
    {
        // Set up to loop over the hits in the first track. 
        Event::McPartToClusPosHitVec trackVec = m_partHitTab->getRelByFirst(mcPart1);
        std::sort(trackVec.begin(),trackVec.end(),CompareTrackHits());

        Event::McPartToClusPosHitVec::const_iterator trackVecIter;
        int                                   bitIndex     = 0;
        int                                   numShared    = 0;

        // Loop over the hits on the track
        for(trackVecIter = trackVec.begin(); trackVecIter != trackVec.end(); trackVecIter++)
        {
            // Recover the cluster associated with this track
            Event::TkrCluster* cluster = ((*trackVecIter)->getSecond())->getFirst();

            // If no cluster skip to the next hit
            if (!cluster) continue;

            // Recover vector of McParticles associated with this cluster
            Event::McPartToClusVec partVec = m_partClusTab->getRelBySecond(cluster);

            // If more than one McParticle associated with the cluster then it is shared
            if (partVec.size() > 1) 
            {
                // Loop through this vector looking for a match to the second track
                for(Event::McPartToClusVec::const_iterator mcPartVecIter = partVec.begin();
                    mcPartVecIter != partVec.end(); mcPartVecIter++)
                {
                    if ((*mcPartVecIter)->getFirst() == mcPart2)
                    {
                        if (numShared < 15) numShared++;
                        else                hitInfo |= 0x08000000;
                        if (bitIndex < 24)  hitInfo |= 1 << bitIndex;
                        else                hitInfo |= 0x00800000;
                    }
                }
            }

            bitIndex++;
        }

        hitInfo = (hitInfo << 4) + numShared;
    }

    return hitInfo;
}

