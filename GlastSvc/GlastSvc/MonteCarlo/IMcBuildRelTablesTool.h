/**
 * @class IMcBuildRelTablesTool
 *
 * @brief Interface to the tool which builds the Monte Carlo event structure, stored in 
 *        the TDS in the McEventStructure object. It also builds the relation table between
 *        McParticles which leave tracks in the Tracker and their associated McPositionHits
 *
 * @author Tracy Usher
 */

#ifndef IMcBuildRelTablesTool_H
#define IMcBuildRelTablesTool_H

#include "GaudiKernel/IAlgTool.h"
#include "Event/MonteCarlo/McEventStructure.h"

static const InterfaceID IID_IMcBuildRelTablesTool("IMcBuildRelTablesTool", 2 , 0);

class IMcBuildRelTablesTool : virtual public IAlgTool 
{
 public:
    /// Retrieve interface ID
    static const InterfaceID& interfaceID() { return IID_IMcBuildRelTablesTool; }

    /// Builds the event structure
    virtual void buildEventStructure()=0;

    /// Builds tables relating Tracker "hits" to McParticles
    virtual void buildMonteCarloTracks()=0;

    /// Builds tables relating MC tracks to TKR Pattern Candidate tracks
    virtual void buildMcPatCandRelations()=0;
};
#endif
