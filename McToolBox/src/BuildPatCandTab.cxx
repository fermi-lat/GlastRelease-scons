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
#include "Event/Recon/TkrRecon/TkrClusterCol.h"
#include "Event/Recon/TkrRecon/TkrPatCand.h"

Event::BuildPatCandTab::BuildPatCandTab(IDataProviderSvc* dataSvc)
{
    StatusCode sc = StatusCode::SUCCESS;
 
    // Look up the Pattern Track collection from the TDS
    Event::TkrPatCandCol* pTkrCands = SmartDataPtr<Event::TkrPatCandCol>(dataSvc,EventModel::TkrRecon::TkrPatCandCol);
    Event::TkrClusterCol* pTkrClus  = SmartDataPtr<Event::TkrClusterCol>(dataSvc,EventModel::TkrRecon::TkrClusterCol); 

    // Look up the Monte Carlo track tables 
    SmartDataPtr<Event::McPartToClusTabList> mcClusTable(dataSvc,EventModel::MC::McPartToClusTab);
    Event::McPartToClusTab mcPartToClusTab(mcClusTable);

    //Create the pattern track relational tables
    Event::McPartToTkrCandHitTab partToHitTab;
    Event::McPartToTkrPatCandTab partToPatCandTab;

    partToHitTab.init();
    partToPatCandTab.init();

    //Store them in the TDS
    sc = dataSvc->registerObject(EventModel::MC::McPartToTkrCandHitTab,partToHitTab.getAllRelations());
    sc = dataSvc->registerObject(EventModel::MC::McPartToTkrPatCandTab,partToPatCandTab.getAllRelations());

    // Loop over the Pattern Candidate tracks
    for(Event::TkrPatCandColPtr cands = pTkrCands->begin(); cands != pTkrCands->end(); cands++)
    {
        Event::TkrPatCand*      patCand = *cands;
        int                     numHits = patCand->numPatCandHits();
        Event::CandHitVectorPtr candPtr = patCand->getHitIterBegin();

        // Loop over the hits in the candidate track
        while(numHits--)
        {
            Event::TkrPatCandHit*  candHit = *candPtr++;
            int                    clusIdx = candHit->HitIndex();
            Event::TkrCluster*     cluster = pTkrClus->getHit(clusIdx);
            Event::McPartToClusVec hitVec  = mcPartToClusTab.getRelBySecond(cluster);

            // Ok, now loop over the  McParticle <-> TkrCluster relations
            Event::McPartToClusVec::const_iterator hitVecIter;
            for(hitVecIter = hitVec.begin(); hitVecIter != hitVec.end(); hitVecIter++)
            {
                Event::McParticle* mcPart = (*hitVecIter)->getFirst();

                Event::McPartToTkrCandHitRel* partToHitRel = new Event::McPartToTkrCandHitRel(mcPart,candHit);
                partToHitRel->addInfo("h");
                if (!partToHitTab.addRelation(partToHitRel)) delete partToHitRel;

                Event::McPartToTkrPatCandRel* partToPatRel = new Event::McPartToTkrPatCandRel(mcPart,patCand);
                partToPatRel->addInfo("p");
                if (!partToPatCandTab.addRelation(partToPatRel)) delete partToPatRel;
            }
        }
    }
    return;    
}

Event::BuildPatCandTab::~BuildPatCandTab()
{
    return;
}
