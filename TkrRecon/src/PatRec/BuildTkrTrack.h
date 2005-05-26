/**
 * @class BuildTkrTrack
 *
 * @brief This class will build a TkrTrack, with TkrTrackHits, given input parameters
 *
 * @author Tracy Usher
 *
 * $Header$
 */

#ifndef BuildTkrTrack_h
#define BuildTkrTrack_h

#include "TkrUtil/ITkrGeometrySvc.h"
#include "Event/Recon/TkrRecon/TkrTrack.h"
#include "src/Track/TkrControl.h"

class BuildTkrTrack 
{
public:
    // Constructors
    BuildTkrTrack(const ITkrGeometrySvc* tkrGeo);

    ~BuildTkrTrack() {}

    Event::TkrTrack* makeNewTkrTrack(Point  startPos, 
                                     Vector startDir, 
                                     double energy, 
                                     std::vector<const Event::TkrCluster*>& clusterVec);

    /// This for adding a single TkrTrackHit given a cluster
    Event::TkrTrackHit* makeTkrTrackHit(const Event::TkrCluster* cluster);

    /// This sets the first hit parameters AFTER track has been constructed
    bool                setFirstHitParams(Event::TkrTrack* track);


private:
    const ITkrGeometrySvc* m_tkrGeom;
    TkrControl*            m_control;
};

#endif

