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
#include "GaudiKernel/ParticleProperty.h"

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
        mcEvent = newMcEventStructure(dataSvc, ppSvc);

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


Event::McEventStructure* Event::BuildEventStructure::newMcEventStructure(IDataProviderSvc* dataSvc, IParticlePropertySvc* partPropSvc)
{
    //Pointer to the class (if it gets created)
    Event::McEventStructure* mcEvent = 0;

    //The primary particle
    Event::McParticleRef primary;

    primary = 0;

    //SmartDataPtr<Event::McParticleCol> mcParts(dataSvc, EventModel::MC::McParticleCol);
    SmartDataPtr<Event::McParticleCol> mcParts(dataSvc, "/Event/MC/McParticleCol");
    Event::McParticleCol::iterator mcPartIter;

    // First step is to find the primary particle
    int numParts = mcParts->size();
    for(mcPartIter = mcParts->begin(); mcPartIter != mcParts->end(); mcPartIter++)
    {
        const Event::McParticle* mcPart = *mcPartIter;

        //The particle we want is a result of the "primary"...
        if (mcPart->getProcess() == std::string("primary"))
        {
            primary = mcPart;
            break;
        }
    }

    // If no primary found then this is an error, but check non zero primary and continue
    if (primary)
    {
        //Attempt to classify the event 
        unsigned long classification = 0;
        int           primaryId      = primary->particleProperty();

        ParticleProperty* ppty     = partPropSvc->findByStdHepID( primaryId );
        std::string       partName = ppty->particle(); 

        // Charged or neutral primary
        if (ppty->charge() == 0)
        {
            classification |= McEventStructure::NEUTRAL;
        }
        else
        {
            classification |= McEventStructure::CHARGED;
        }

        // Label as a gamma as well (if that is what it is...)
        if (partName == "gamma") classification |= McEventStructure::GAMMA;

        // Create the McEventStructure object
        mcEvent = new Event::McEventStructure();

        mcEvent->setPrimaryParticle(primary);

        // Set up to loop over primary daughters 
        const SmartRefVector<Event::McParticle>& daughterVec = primary->daughterList();
        SmartRefVector<Event::McParticle>::const_iterator daughterVecIter;

        // Loop over daughters looking for secondary tracks (hits in the tracker) 
        // Also set process bits if a gamma
        for(daughterVecIter = daughterVec.begin(); daughterVecIter != daughterVec.end(); daughterVecIter++)
        {
            const Event::McParticle* mcPart   = *daughterVecIter;
            bool                     isTrkHit = ((mcPart->statusFlags() & Event::McParticle::POSHIT)) != 0;

            // If particle is a gamma then set the decay process that produced the daughter
            if (partName == "gamma")
            {
                const std::string        process  = mcPart->getProcess();

                if (process == "conv")
                {
                    classification |= McEventStructure::CONVERT;
                    if (isTrkHit) classification |= McEventStructure::TRKCONVERT;
                }
                else if (process == "brem")
                {
                    classification |= McEventStructure::BREM;
                    if (isTrkHit) classification |= McEventStructure::TRKBREM;
                }
                else if (process == "compt")
                {
                    classification |= McEventStructure::COMPT;
                    if (isTrkHit) classification |= McEventStructure::TRKCOMPT;
                }
                else if (process == "phot")
                {
                    classification |= McEventStructure::PHOT;
                    if (isTrkHit) classification |= McEventStructure::TRKPHOT;
                }
                else
                {
                    classification |= McEventStructure::OTHER;
                    if (isTrkHit) classification |= McEventStructure::TRKOTHER;
                }
            }

            // If secondary track add to the list
            if (isTrkHit)
            {
                // Add the secondary to the list
                mcEvent->addSecondary(mcPart);

                // Go through and find all the tracks "associated" to this one
                findAssociated(mcEvent, mcPart);
            }
        }

        // Update the classification bits
        mcEvent->setClassificationBits(classification);
    }

    int numScndrys = mcEvent->getNumSecondaries();
    int numAssoc   = mcEvent->getNumAssociated();
  
    return mcEvent;    
}
/*
bool Event::BuildEventStructure::isPrimaryDaughter(const Event::McParticle* primary, const Event::McParticle* mcPart)
{
    //Search the primary particles daughter list for a match
    const SmartRefVector<Event::McParticle>& daughterVec = primary->daughterList();
    SmartRefVector<Event::McParticle>::const_iterator daughterVecIter;

    for(daughterVecIter = daughterVec.begin(); daughterVecIter != daughterVec.end(); daughterVecIter++)
    {
        //If a match then return true
        if (mcPart == *daughterVecIter) return true;
    }

    return false;
}
*/
void Event::BuildEventStructure::findAssociated(Event::McEventStructure* mcEvent, const Event::McParticle* mcPart)
{
    //Search the current particle's daughter list for particles which leave McPositionHits
    const SmartRefVector<Event::McParticle>& daughterVec = mcPart->daughterList();
    SmartRefVector<Event::McParticle>::const_iterator daughterVecIter;

    for(daughterVecIter = daughterVec.begin(); daughterVecIter != daughterVec.end(); daughterVecIter++)
    {
        const Event::McParticle* daughter = *daughterVecIter;

        // If it leaves an McPositionHit then add to the list
        if (daughter->statusFlags() & Event::McParticle::POSHIT) mcEvent->addAssociated(daughter);

        // Look at this particle's daughters to see if they leave hits...
        idents::VolumeIdentifier iniVolId = const_cast<Event::McParticle*>(mcPart)->getInitialId();
        idents::VolumeIdentifier finVolId = const_cast<Event::McParticle*>(mcPart)->getFinalId();

        // If the particle starts or stops in the tracker, then look at its daughters
        if ((iniVolId[0] == 0 && iniVolId[3] == 1 && iniVolId.size() > 3) ||
            (finVolId[0] == 0 && finVolId[3] == 1 && finVolId.size() > 3)   )
        {
            findAssociated(mcEvent, daughter);
        }
    }

    return;
}
