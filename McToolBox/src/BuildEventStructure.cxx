 /**
 * @class BuildEventStructure
 *
 * @brief Implements a Gaudi Tool for setting the candidate track energies before 
 *        the track fit
 *
 * @author The Tracking Software Group
 *
 * $Header$
 */
#include "BuildEventStructure.h"
#include "Event/MonteCarlo/McRelTableDefs.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/GaudiException.h" 
#include "Event/TopLevel/EventModel.h"

Event::BuildEventStructure::BuildEventStructure(IDataProviderSvc* dataSvc, IParticlePropertySvc* ppSvc)
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

    // Retrieve the pointer to the McEventStructure
    SmartDataPtr<Event::McEventStructure> mcEvent(dataSvc,EventModel::MC::McEventStructure);

    // If it doesn't exist then we need to build the MC structure
    if (mcEvent == 0)
    {
        //This builds the Monte Carlo event structure - basically a description of the event
        mcEvent = new Event::McEventStructure(dataSvc, ppSvc);

        //Store in the Mci TDS 
        sc = dataSvc->registerObject(EventModel::MC::McEventStructure, mcEvent);
        if (sc.isFailure())
        {
            throw GaudiException("Cannot store the McEventStructure in the TDS", "BuildEventStructure", sc);
        }

        //Define and initialize the table relating an McParticle to its McPositionHits
        Event::McPartToPosHitTab  mcPartToPosHitTab;
        mcPartToPosHitTab.init();

        //Store the table in the TDS
        sc = dataSvc->registerObject(EventModel::MC::McPartToPosHitTab,mcPartToPosHitTab.getAllRelations());
        if (sc.isFailure())
        {
            throw GaudiException("Cannot store the McPartToPosHitTab in the TDS", "BuildEventStructure", sc);
        }

        //Retrieve the collection of McPositionHits from the TDS
        SmartDataPtr<Event::McPositionHitVector> mcPosHits(dataSvc, EventModel::MC::McPositionHitCol);
        Event::McPositionHitVector::iterator mcPosHitIter;

        //Get the vector of "valid" McParticles
        Event::McParticleRefVec mcPartVec = mcEvent->getTrackVector();

        //Loop through McPositionHits and make relations from "valid" McParticle matches
        for(mcPosHitIter = mcPosHits->begin(); mcPosHitIter != mcPosHits->end(); mcPosHitIter++)
        {
            const Event::McPositionHit* mcPosHit = *mcPosHitIter;
            const Event::McParticle*    mcPart   = mcPosHit->mcParticle();

            // Make sure this McPositionHit is related to one of our McEventStructure particles
            if (find(mcPartVec.begin(),mcPartVec.end(),mcPart))
            {
                Event::McPartToPosHitRel* posHitRel = 
                    new Event::McPartToPosHitRel(const_cast<Event::McParticle*>(mcPart), const_cast<Event::McPositionHit*>(mcPosHit));
                
                if (!mcPartToPosHitTab.addRelation(posHitRel)) delete posHitRel;
            }
        }

    }

    // finished
    return;
}

bool Event::BuildEventStructure::find(const Event::McParticleRefVec::iterator begin, 
                                      const Event::McParticleRefVec::iterator end,
                                      const Event::McParticle*          mcPart)
{
    for(Event::McParticleRefVec::iterator iter = begin; iter != end; iter++)
    {
        const Event::McParticle* testMcPart = *iter;

        if (testMcPart == mcPart) return true;
    }

    return false;
}


Event::BuildEventStructure::~BuildEventStructure()
{
    return;
}
