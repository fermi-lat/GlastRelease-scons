/**
 * @class McEventStructure
 *
 * @brief This object attempts to use the Monte Carlo information to organize and
 *        classify a single event in LAT. To that end, it identifies the primary 
 *        particle and its interactions within the LAT. After this it identifies any
 *        direct daughters of the primary particle WHICH RESULT IN MCPOSITIONHITs in
 *        the tracker, and then identifies any "associated" (daughters of daughters) 
 *        which also leave McPositionHits in the tracker. 
 *
 *        The type of event, charged or gamma, the type of interaction in the LAT (also
 *        within the tracker) can be determined from the event classification bit string.
 *
 * @author Tracy Usher
 *
 * $Header$
 */

#ifndef McEventStructure_h
#define McEventStructure_h

#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/IParticlePropertySvc.h"
#include "GaudiKernel/ContainedObject.h"
#include "GaudiKernel/SmartRefVector.h"
#include "Event/MonteCarlo/McParticle.h"
#include "Event/MonteCarlo/McPositionHit.h"
#include "Event/RelTable/RelTable.h"

namespace Event {

typedef SmartRefVector<Event::McParticle> McParticleRefVec; 
typedef SmartRef<Event::McParticle>       McParticleRef;

class McEventStructure : public DataObject
{
public:
    //! Define bits to help classify the event
    enum ClassificationBits{  
        NOPRIMARY  = 1 ,    //! No primary particle found (can't happen?)
        CHARGED    = 1<<1,  //! Primary particle is charged
        NEUTRAL    = 1<<2,  //! Primary particle is neutral
        GAMMA      = 1<<3,  //! Primary is a gamma
        CONVERT    = 1<<4,  //! Secondaries from gamma conversion
        BREM       = 1<<5,  //! Secondaries from Bremstrahlung
        COMPT      = 1<<6,  //! Secondaries from Compton scatter
        PHOT       = 1<<7,  //! Secondaries from Photoelectric effect (?)
        OTHER      = 1<<8,  //! Secondaries from some other process 
        TRKCONVERT = 1<<12, //! Secondaries from gamma conversion in tracker
        TRKBREM    = 1<<13, //! Secondaries from Bremstrahlung in tracker
        TRKCOMPT   = 1<<14, //! Secondaries from Compton scatter in tracker
        TRKPHOT    = 1<<15, //! Secondaries from Photoelectric effect in tracker
        TRKOTHER   = 1<<16, //! Secondaries from some other process in tracker
        RUNBIT     = 1<<24  //! Indicates code was called (does nothing)
    };

    //! Dataobject compliant constructor
    McEventStructure(IDataProviderSvc* dataSvc, IParticlePropertySvc* partPropSvc);
   ~McEventStructure();

    //! Retrieve classification bits (see above definitions)
    const unsigned long                     getClassificationBits() const {return m_classification;}
    
    //! Retrieve reference to the primary particle
    const Event::McParticleRef              getPrimaryParticle()    const {return m_primary;}

    //! Retrieve number and iterators to the daughters of the primary with hits in tracker
    const int                               getNumSecondaries()     const {return m_secondaries.size();}
    Event::McParticleRefVec::const_iterator beginSecondaries()      const {return m_secondaries.begin();}
    Event::McParticleRefVec::const_iterator endSecondaries()        const {return m_secondaries.end();}

    //! Retrieve number and iterators to the particles associated with primary 
    const int                               getNumAssociated()      const {return m_associated.size();}
    Event::McParticleRefVec::const_iterator beginAssociated()       const {return m_associated.begin();}
    Event::McParticleRefVec::const_iterator endAssociated()         const {return m_associated.end();}

    //! Return an McParticle reference vector of all McParticles which leave hits in tracker
    Event::McParticleRefVec                 getTrackVector();

private:
    bool isPrimaryDaughter(const Event::McParticle* mcPart);
    void findAssociated(const Event::McParticle* mcPart);

    /// Bit-field for classification
    unsigned long             m_classification;

    Event::McParticleRef      m_primary;
    Event::McParticleRefVec   m_secondaries;
    Event::McParticleRefVec   m_associated;
};

/*
// typedefs for relating McParticles (above) to associated McPositionHits
typedef Event::RelTable<Event::McParticle, Event::McPositionHit>   McPartToPosHitTab;
typedef Event::Relation<Event::McParticle, Event::McPositionHit>   McPartToPosHitRel;
typedef ObjectList<McPartToPosHitRel>                              McPartToPosHitTabList;
typedef std::vector<Event::McPartToPosHitRel*>                     McPartToPosHitVec;
*/
};

#endif