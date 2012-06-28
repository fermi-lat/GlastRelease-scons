/** @file ITkrTreeTrackFinder.h
 * @class ITkrTreeTrackFinder
 *
 * @brief This provides an interface to access the code which extracts tracks from trees
 *
 * @author Tracy Usher
 *
 * $Header$
 *
*/

#ifndef __ITkrTreeTrackFinder_H
#define __ITkrTreeTrackFinder_H 1

#include "GaudiKernel/IAlgTool.h"
#include "GaudiKernel/IProperty.h"

#include "Event/Recon/TkrRecon/TkrTree.h"
#include "Event/Recon/CalRecon/CalCluster.h"

static const InterfaceID IID_ITkrTreeTrackFinder("ITkrTreeTrackFinder", 7111 , 0);

class ITkrTreeTrackFinder : virtual public IAlgTool
{
public:
    /// Retrieve interface ID
    static const InterfaceID& interfaceID() { return IID_ITkrTreeTrackFinder; }

    /// Method called to find the tracks for a given tree
    virtual int findTracks(Event::TkrTree* tree, double eventEnergy, Event::CalCluster* cluster = 0) = 0;
};

#endif
