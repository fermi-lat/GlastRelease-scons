/**
 * @class IMcGetTrackInfoTool
 *
 * @brief Interface to the tool for returning information from the tables relating Monte Carlo 
 *        McParticles to TkrPatCand and TkrPatCandHit objects. This information can be 
 *        used to study the performance of the pattern recognition.
 *
 * @author Tracy Usher
 */

#ifndef IMcGetTrackInfoTool_H
#define IMcGetTrackInfoTool_H

#include "GaudiKernel/IAlgTool.h"
#include "Event/MonteCarlo/McEventStructure.h"
#include "Event/MonteCarlo/McRelTableDefs.h"
#include "Event/Recon/TkrRecon/TkrPatCand.h"

static const InterfaceID IID_IMcGetTrackInfoTool("IMcGetTrackInfoTool", 2 , 0);

class IMcGetTrackInfoTool : virtual public IAlgTool 
{
 public:
    /// Retrieve interface ID
    static const InterfaceID& interfaceID() { return IID_IMcGetTrackInfoTool; }

    /// @brief Return the number of McParticles associated to a given TkrPatCand track
    virtual int                           getNMcParticles(const Event::TkrPatCand*)=0;

    /// @brief Return the vector of McParticles associated to a given TkrPatCand track
    virtual const Event::McParticleRefVec getMcParticleVec(const Event::TkrPatCand*)=0;

    /// @brief Return TkrPatCand <-> McParticle "best" match
    virtual const Event::McParticle*      getBestMcParticle(const Event::TkrPatCand*)=0;

    /// @brief Return TkrPatCand <-> McParticle "best" match
    virtual const Event::TkrPatCand*      getBestTkrPatCand(const Event::McParticle*)=0;

    /// @brief Return the number of hits on a given TkrPatCand track associated to
    ///        a given McParticle
    virtual int                           getNumMcHits(const Event::TkrPatCand*, const Event::McParticle*)=0;
};
#endif
