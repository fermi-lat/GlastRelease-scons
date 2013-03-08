 /**
 * @class BuildMcTracks
 *
 * @brief This object builds tables to relate McParticles, McPositionHits and TkrClusters
 *        together in order to associate hits in the Tracker with Monte Carlo tracks.
 *        Effectively, this association can be made with four relational tables:
 *        1) A table relating McParticles with McPositionHits    (built in BuildEventStructure)
 *        2) A table relating McPositionHits with TkrClusters    (built in TkrRecon)
 *        3) A table relating McParticles with TkrClusters
 *        4) The most complicated, a table relating McParticles with existing 
 *           TkrCluster<->McPositionHit relations
 *
 *        So, here we only need build the last two tables (piece of cake!)
 *
 * @author The Tracking Software Group
 *
 * $Header$
 */
#include "BuildMcTracks.h"
#include "Event/MonteCarlo/McRelTableDefs.h"
#include "Event/MonteCarlo/McParticle.h"
#include "Event/MonteCarlo/McPositionHit.h"
#include "Event/Recon/TkrRecon/TkrCluster.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "Event/TopLevel/EventModel.h"
#include "idents/TkrId.h"
#include <stdexcept>

Event::BuildMcTracks::BuildMcTracks(IDataProviderSvc* dataSvc)
{
    //Take care of insuring that data area has been created
    StatusCode  sc;

    //Define the table needed for the MC tracks
    Event::McPartToClusTab       mcPartToClusTab;
    Event::McPartToClusPosHitTab mcPartToClusPosHitTab;

    //Initialize the tables
    mcPartToClusTab.init();
    mcPartToClusPosHitTab.init();

    //Store these objects in the TDS
    sc = dataSvc->registerObject(EventModel::MC::McPartToClusTab,mcPartToClusTab.getAllRelations());
    sc = dataSvc->registerObject(EventModel::MC::McPartToClusHitTab,mcPartToClusPosHitTab.getAllRelations());

    //Recover the McPositionHit to Cluster relational table
    SmartDataPtr<Event::ClusMcPosHitTabList> tkrTable(dataSvc,EventModel::Digi::TkrClusterHitTab);
    Event::ClusMcPosHitTab clusHitTab(tkrTable);

    //If no relations then no reason to continue
    if ((clusHitTab.getAllRelations())->size() == 0) return;
 
    // Get the MC info we want from the TDS
    SmartDataPtr<Event::McPositionHitVector> posHits(dataSvc, EventModel::MC::McPositionHitCol);

    // Loop through McPositionHits to build McSiLayerHits 
    Event::McPositionHitVector::const_iterator hit;
    for (hit = posHits->begin(); hit != posHits->end(); hit++ ) 
    {
        const Event::McPositionHit*    mcPosHit = *hit;
        const Event::McParticle*       mcPart   = mcPosHit->mcParticle();
        const idents::VolumeIdentifier curVolId = mcPosHit->volumeID();

        // Want only hits directly from parent
        if (mcPosHit->getPackedFlags() & 0x10000000)
        {
            continue;
        }

        try
        {
            idents::TkrId trkId(curVolId);
        }
        catch (std::invalid_argument ex) 
        {
            continue;
        }

        // Only interested in volumes inside the tracker
        if (curVolId[0] != 0) continue;

        // Find the cluster associated with this McPositionHit
        std::vector<Event::ClusMcPosHitRel*> clusRelVector = clusHitTab.getRelBySecond(mcPosHit);

        // If no clusters associated to this McPositionHit, then skip
        int numRels = clusRelVector.size();

        if (numRels == 0) continue;

        // Get the cluster associated with this McPositionHit
        Event::ClusMcPosHitRel* clusMcPosHitRel = numRels > 0 ? clusRelVector[0] : 0;
        Event::TkrCluster*      tkrCluster      = numRels > 0 ? clusRelVector[0]->getFirst() : 0;

        // Patch up the clusHitTab (temporarily done here)
        if (numRels == 0)
        {
            clusMcPosHitRel = new Event::ClusMcPosHitRel(tkrCluster, const_cast<Event::McPositionHit*>(mcPosHit));
            if (!clusHitTab.addRelation(clusMcPosHitRel)) delete clusMcPosHitRel;
        }

        // Create the McParticle to TkrCluster relation and add to the table
        Event::McPartToClusRel* partClusRel = new Event::McPartToClusRel(const_cast<Event::McParticle*>(mcPart), tkrCluster);
        if (!mcPartToClusTab.addRelation(partClusRel)) delete partClusRel;

        Event::McPartToClusPosHitRel* partClusHitRel = new Event::McPartToClusPosHitRel(const_cast<Event::McParticle*>(mcPart), clusMcPosHitRel);
        if (!mcPartToClusPosHitTab.addRelation(partClusHitRel)) delete partClusHitRel;
    }

    return;
}

Event::BuildMcTracks::~BuildMcTracks()
{
    return;
}
