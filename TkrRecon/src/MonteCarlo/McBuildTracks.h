/**
 * @class McBuildTracks
 *
 * @brief Represents a Monte Carlo hit in a tracker layer
 *
 * @author The Tracking Software Group
 *
 * $Header$
 */
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/ContainedObject.h"
#include "GaudiKernel/SmartRefVector.h"
#include "Event/MonteCarlo/McParticle.h"

#ifndef McBuildTracks_h
#define McBuildTracks_h

namespace Event {

class McBuildTracks : virtual public ContainedObject
{
public:
    /// Standard Gaudi Tool interface constructor
    McBuildTracks(IDataProviderSvc* dataSvc);
   ~McBuildTracks();

private:
};

};

#endif