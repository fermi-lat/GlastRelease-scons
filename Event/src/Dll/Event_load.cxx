/** 
* @file TkrUtil_load.cxx
* @brief This is needed for forcing the linker to load all components
* of the library.
*
*  $Header$
*/

#include "GaudiKernel/DeclareFactoryEntries.h"

// Define a hook for declaring Contained objects
// Why doesn't Gaudi define this?
#define DECLARE_CONTAINEDOBJECT(x) extern const IContainedObjectFactory& x##Factory; x##Factory.clID();

// This will get around our namespace problems
#define DECLARE_NAMESPACE_DATAOBJECT(n,x) extern const IDataObjectFactory& n##x##Factory; n##x##Factory.clID();

#define DECLARE_NAMESPACE_DATAOBJECT_FACTORY( n,x )     \
static DataObjectFactory< n::x > s_##n##x##Factory; \
const IDataObjectFactory& n##x##Factory = s_##n##x##Factory;

//
// This declares the top level objects in the TDS
//
#include "Event/TopLevel/Event.h"
#include "Event/TopLevel/MCEvent.h"
#include "Event/TopLevel/DigiEvent.h"

using namespace Event; // Makes life easier...

DECLARE_DATAOBJECT_FACTORY(EventHeader);
DECLARE_DATAOBJECT_FACTORY(MCEvent);
DECLARE_DATAOBJECT_FACTORY(DigiEvent);

//
// This declares the Monte Carlo TDS objects
//
#include "Event/MonteCarlo/D2Entry.h"
#include "Event/MonteCarlo/Exposure.h"
#include "Event/MonteCarlo/McIntegratingHit.h"
#include "Event/MonteCarlo/McParticle.h"
#include "Event/MonteCarlo/McPositionHit.h"
#include "Event/MonteCarlo/McTkrStrip.h"
#include "Event/MonteCarlo/McTrajectory.h"
DECLARE_CONTAINEDOBJECT_FACTORY(D2Entry);
DECLARE_CONTAINEDOBJECT_FACTORY(Exposure);
DECLARE_CONTAINEDOBJECT_FACTORY(McIntegratingHit);
DECLARE_CONTAINEDOBJECT_FACTORY(McParticle);
DECLARE_CONTAINEDOBJECT_FACTORY(McPositionHit);
DECLARE_CONTAINEDOBJECT_FACTORY(McTkrStrip);
DECLARE_CONTAINEDOBJECT_FACTORY(McTrajectory);

//
// This declares the Digi TDS objects
//
#include "Event/Digi/AcdDigi.h"
#include "Event/Digi/CalDigi.h"
#include "Event/Digi/GltDigi.h"
#include "Event/Digi/TkrDigi.h"
DECLARE_CONTAINEDOBJECT_FACTORY(AcdDigi);
DECLARE_CONTAINEDOBJECT_FACTORY(CalDigi);
DECLARE_CONTAINEDOBJECT_FACTORY(GltDigi);
DECLARE_CONTAINEDOBJECT_FACTORY(TkrDigi);

//
// This declares the Recon TDS objects
#include "Event/Recon/AcdRecon/AcdRecon.h"
#include "Event/Recon/CalRecon/CalRecon.h"
#include "Event/Recon/CalRecon/CalXtalRecData.h"
DECLARE_DATAOBJECT_FACTORY(AcdRecon);
DECLARE_DATAOBJECT_FACTORY(CalRecon);
DECLARE_CONTAINEDOBJECT_FACTORY(CalXtalRecData);

//
// Do tracker stuff separately
//
#include "Event/Recon/TkrRecon/TkrCluster.h"
#include "Event/Recon/TkrRecon/TkrFitTrack.h"
#include "Event/Recon/TkrRecon/TkrFitTrackBase.h"
#include "Event/Recon/TkrRecon/TkrPatCand.h"
#include "Event/Recon/TkrRecon/TkrPatCandHit.h"
#include "Event/Recon/TkrRecon/TkrRecon.h"
#include "Event/Recon/TkrRecon/TkrTrack.h"
#include "Event/Recon/TkrRecon/TkrTrackHit.h"
#include "Event/Recon/TkrRecon/TkrVertex.h"
DECLARE_CONTAINEDOBJECT_FACTORY(TkrCluster);
DECLARE_CONTAINEDOBJECT_FACTORY(TkrFitTrack);
DECLARE_CONTAINEDOBJECT_FACTORY(TkrPatCand);
DECLARE_CONTAINEDOBJECT_FACTORY(TkrPatCandHit);
DECLARE_CONTAINEDOBJECT_FACTORY(TkrTrack);
DECLARE_CONTAINEDOBJECT_FACTORY(TkrTrackHit);
DECLARE_CONTAINEDOBJECT_FACTORY(TkrVertex);
DECLARE_DATAOBJECT_FACTORY(TkrClusterCol);
DECLARE_DATAOBJECT_FACTORY(TkrRecon);

DECLARE_FACTORY_ENTRIES(Event) 
{
    DECLARE_DATAOBJECT( EventHeader  );
    DECLARE_DATAOBJECT( MCEvent      );
    DECLARE_DATAOBJECT( DigiEvent    );

    DECLARE_CONTAINEDOBJECT( D2Entry          );
    DECLARE_CONTAINEDOBJECT( Exposure         );
    DECLARE_CONTAINEDOBJECT( McIntegratingHit );
    DECLARE_CONTAINEDOBJECT( McParticle       );
    DECLARE_CONTAINEDOBJECT( McPositionHit    );
    DECLARE_CONTAINEDOBJECT( McTkrStrip       );
    DECLARE_CONTAINEDOBJECT( McTrajectory     );

    DECLARE_CONTAINEDOBJECT( AcdDigi          );
    DECLARE_CONTAINEDOBJECT( CalDigi          );
    DECLARE_CONTAINEDOBJECT( GltDigi          );
    DECLARE_CONTAINEDOBJECT( TkrDigi          );

    DECLARE_DATAOBJECT(      AcdRecon         );
    DECLARE_DATAOBJECT(      CalRecon         );
    DECLARE_CONTAINEDOBJECT( CalXtalRecData   );

    DECLARE_CONTAINEDOBJECT( TkrCluster       );
    DECLARE_CONTAINEDOBJECT( TkrFitTrack      );
    DECLARE_CONTAINEDOBJECT( TkrPatCand       );
    DECLARE_CONTAINEDOBJECT( TkrPatCandHit    );
    DECLARE_CONTAINEDOBJECT( TkrTrack         );
    DECLARE_CONTAINEDOBJECT( TkrTrackHit      );
    DECLARE_CONTAINEDOBJECT( TkrVertex        );

    DECLARE_DATAOBJECT(      TkrRecon         );
    DECLARE_DATAOBJECT(      TkrClusterCol    );
} 


