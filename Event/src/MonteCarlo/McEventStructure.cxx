/**
 * @class McEventStructure
 *
 * @brief Implements a Gaudi Tool for setting the candidate track energies before 
 *        the track fit
 *
 * @author The Tracking Software Group
 *
 * $Header$
 */
#include "Event/MonteCarlo/McEventStructure.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/ParticleProperty.h"
#include "Event/TopLevel/EventModel.h"


Event::McEventStructure::McEventStructure(IDataProviderSvc* dataSvc, IParticlePropertySvc* partPropSvc) :
                                          m_classification(McEventStructure::RUNBIT)
{
    // clear refvec's (just in case)
    m_primary = 0;
    m_secondaries.clear();
    m_associated.clear();

    SmartDataPtr<Event::McParticleCol> mcParts(dataSvc, EventModel::MC::McParticleCol);
    Event::McParticleCol::iterator mcPartIter;

    // First step is to find the primary particle
    int numParts = mcParts->size();
    for(mcPartIter = mcParts->begin(); mcPartIter != mcParts->end(); mcPartIter++)
    {
        const Event::McParticle* mcPart = *mcPartIter;

        //The particle we want is a result of the "primary"...
        if (mcPart->getProcess() == std::string("primary"))
        {
            m_primary = mcPart;
            break;
        }
    }

    // If no primary found then this is an error, but check non zero primary and continue
    if (m_primary)
    {
        //Attempt to classify the event 
        int primaryId = m_primary->particleProperty();

        ParticleProperty* ppty     = partPropSvc->findByStdHepID( primaryId );
        std::string       partName = ppty->particle(); 

        // Charged or neutral primary
        if (ppty->charge() == 0)
        {
            m_classification |= McEventStructure::NEUTRAL;
        }
        else
        {
            m_classification |= McEventStructure::CHARGED;
        }

        // Label as a gamma as well (if that is what it is...)
        if (partName == "gamma") m_classification |= McEventStructure::GAMMA;

        // Set up to loop over primary daughters 
        const SmartRefVector<Event::McParticle>& daughterVec = m_primary->daughterList();
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
                    m_classification |= McEventStructure::CONVERT;
                    if (isTrkHit) m_classification |= McEventStructure::TRKCONVERT;
                }
                else if (process == "brem")
                {
                    m_classification |= McEventStructure::BREM;
                    if (isTrkHit) m_classification |= McEventStructure::TRKBREM;
                }
                else if (process == "compt")
                {
                    m_classification |= McEventStructure::COMPT;
                    if (isTrkHit) m_classification |= McEventStructure::TRKCOMPT;
                }
                else if (process == "phot")
                {
                    m_classification |= McEventStructure::PHOT;
                    if (isTrkHit) m_classification |= McEventStructure::TRKPHOT;
                }
                else
                {
                    m_classification |= McEventStructure::OTHER;
                    if (isTrkHit) m_classification |= McEventStructure::TRKOTHER;
                }
            }

            // If secondary track add to the list
            if (isTrkHit)
            {
                // Add the secondary to the list
                m_secondaries.push_back(mcPart);

                // Go through and find all the tracks "associated" to this one
                findAssociated(mcPart);
            }
        }

    }

    int numScndrys = m_secondaries.size();
    int numAssoc   = m_associated.size();
  
    return;    
}

Event::McEventStructure::~McEventStructure()
{
    return;
}

bool Event::McEventStructure::isPrimaryDaughter(const Event::McParticle* mcPart)
{
    //Search the primary particles daughter list for a match
    const SmartRefVector<Event::McParticle>& daughterVec = m_primary->daughterList();
    SmartRefVector<Event::McParticle>::const_iterator daughterVecIter;

    for(daughterVecIter = daughterVec.begin(); daughterVecIter != daughterVec.end(); daughterVecIter++)
    {
        //If a match then return true
        if (mcPart == *daughterVecIter) return true;
    }

    return false;
}

void Event::McEventStructure::findAssociated(const Event::McParticle* mcPart)
{
    //Search the current particle's daughter list for particles which leave McPositionHits
    const SmartRefVector<Event::McParticle>& daughterVec = mcPart->daughterList();
    SmartRefVector<Event::McParticle>::const_iterator daughterVecIter;

    for(daughterVecIter = daughterVec.begin(); daughterVecIter != daughterVec.end(); daughterVecIter++)
    {
        const Event::McParticle* daughter = *daughterVecIter;

        // If it leaves an McPositionHit then add to the list
        if (daughter->statusFlags() & Event::McParticle::POSHIT) m_associated.push_back(daughter);

        // Look at this particle's daughters to see if they leave hits...
        idents::VolumeIdentifier iniVolId = const_cast<Event::McParticle*>(mcPart)->getInitialId();
        idents::VolumeIdentifier finVolId = const_cast<Event::McParticle*>(mcPart)->getFinalId();

        // If the particle starts or stops in the tracker, then look at its daughters
        if ((iniVolId[0] == 0 && iniVolId[3] == 1 && iniVolId.size() > 3) ||
            (finVolId[0] == 0 && finVolId[3] == 1 && finVolId.size() > 3)   )
        {
            findAssociated(daughter);
        }
    }

    return;
}

Event::McParticleRefVec Event::McEventStructure::getTrackVector()
{
    Event::McParticleRefVec trackVec;

    trackVec.clear();

    if (m_primary->statusFlags() & Event::McParticle::POSHIT) trackVec.push_back(m_primary);

    Event::McParticleRefVec::const_iterator refIter;

    for(refIter = m_secondaries.begin(); refIter != m_secondaries.end(); refIter++) trackVec.push_back(*refIter);
    for(refIter = m_associated.begin();  refIter != m_associated.end();  refIter++) trackVec.push_back(*refIter);

    return trackVec;
}


