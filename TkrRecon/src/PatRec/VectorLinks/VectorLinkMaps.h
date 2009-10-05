/**
 * @class VectorLinkMaps
 *
 * @brief Implementation of a particular "matrix" class for the generic Kalman Filter fitter. 
 *        This version based on CLHEP HepMatrix. 
 *
 * @author Tracy Usher
 *
 * $Header$
 */

#ifndef VectorLinkMaps_h
#define VectorLinkMaps_h

#include "Event/Recon/TkrRecon/TkrTrackElements.h"
#include "Event/Recon/TkrRecon/TkrVecPointsLink.h"
#include "Event/Recon/TkrRecon/TkrCluster.h"

/// Define a map between a TrackElement and all associated TkrPointsLinks 
/// In this case, TrackElement is the key in this map
typedef std::map< const Event::TkrTrackElements*, Event::TkrVecPointsLinkPtrVec> TrackElementToLinksMap;
typedef std::pair<const Event::TkrTrackElements*, Event::TkrVecPointsLinkPtrVec> TrackElementToLinkPair;

/// Define a reverse map taking us from a TrkPointsLink to all TrackElements it 
/// is associated with. 
/// In this case, TkrPointsLink is the key in this map
typedef std::map< const Event::TkrVecPointsLink*, Event::TkrTrackElementsPtrVec> VecLinksToElementsMap;
typedef std::pair<const Event::TkrVecPointsLink*, Event::TkrTrackElementsPtrVec> VecLinksToelementsPair;

#endif

