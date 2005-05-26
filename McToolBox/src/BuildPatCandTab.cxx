/**
 * @class BuildPatCandTab
 *
 * @brief This object builds the relational tables associating the Monte Carlo with the
 *        TkrRecon pattern recognition track candidates. Basically, two tables are created:
 *        1) A table which relates McParticles to TkrPatCands. Because a TkrPatCand can pick
 *           up the "wrong" hits, there are multiple associations between both
 *        2) A Table which relates McParticles to TkrPatCandHits. Since a cluster can be 
 *           formed from several McPositionHits, it is possible for a given track hit to
 *           be associated with more than one McParticle. 
 *
 * @author The Tracking Software Group
 *
 * $Header$
 */
#include "BuildPatCandTab.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "Event/TopLevel/EventModel.h"
#include "Event/MonteCarlo/McRelTableDefs.h"
#include "Event/MonteCarlo/McEventStructure.h"
#include "Event/Recon/TkrRecon/TkrCluster.h"
#include "Event/Recon/TkrRecon/TkrTrack.h"

Event::BuildPatCandTab::BuildPatCandTab(IDataProviderSvc* dataSvc)
{
    StatusCode sc = StatusCode::SUCCESS;
 
    // Look up the Pattern Track collection from the TDS
    Event::TkrTrackCol* pTkrCands = SmartDataPtr<Event::TkrTrackCol>(dataSvc,EventModel::TkrRecon::TkrTrackCol);
    //Event::TkrClusterCol* pTkrClus  = SmartDataPtr<Event::TkrClusterCol>(dataSvc,EventModel::TkrRecon::TkrClusterCol); 

    // Look up the Monte Carlo track tables 
    SmartDataPtr<Event::McPartToClusTabList> mcClusTable(dataSvc,EventModel::MC::McPartToClusTab);
    Event::McPartToClusTab mcPartToClusTab(mcClusTable);

    //Create the pattern track relational tables
    Event::McPartToTkrTrackHitTab partToHitTab;
    Event::McPartToTkrTrackTab    partToTkrTrackTab;

    partToHitTab.init();
    partToTkrTrackTab.init();

    //Store them in the TDS
    sc = dataSvc->registerObject(EventModel::MC::McPartToTkrTrackHitTab,partToHitTab.getAllRelations());
    sc = dataSvc->registerObject(EventModel::MC::McPartToTkrTrackTab,partToTkrTrackTab.getAllRelations());

    // Loop over the Pattern Candidate tracks
    for(Event::TkrTrackColPtr cands = pTkrCands->begin(); cands != pTkrCands->end(); cands++)
    {
        Event::TkrTrack* patCand = *cands;

        // Loop over the hits in the candidate track
        for(Event::TkrTrackHitVecItr hitIter = patCand->begin(); hitIter != patCand->end(); hitIter++)
        {
            Event::TkrTrackHit*      candHit = *hitIter;
            const Event::TkrCluster* cluster = candHit->getClusterPtr();

            if (cluster == 0) continue;

            Event::McPartToClusVec   hitVec  = mcPartToClusTab.getRelBySecond(cluster);

            // Ok, now loop over the  McParticle <-> TkrCluster relations
            Event::McPartToClusVec::const_iterator hitVecIter;
            for(hitVecIter = hitVec.begin(); hitVecIter != hitVec.end(); hitVecIter++)
            {
                Event::McParticle* mcPart = (*hitVecIter)->getFirst();

                Event::McPartToTkrTrackHitRel* partToHitRel = new Event::McPartToTkrTrackHitRel(mcPart,candHit);
                partToHitRel->addInfo("h");
                if (!partToHitTab.addRelation(partToHitRel)) delete partToHitRel;

                Event::McPartToTkrTrackRel* partToTrackRel = new Event::McPartToTkrTrackRel(mcPart,patCand);
                partToTrackRel->addInfo("p");
                if (!partToTkrTrackTab.addRelation(partToTrackRel)) delete partToTrackRel;
            }
        }
    }
    return;    
}

Event::BuildPatCandTab::~BuildPatCandTab()
{
    return;
}
