/** @file TkrTreeTrackFinder.h
 * @class TkrTreeTrackFinder
 *
 * @brief Provides X-Y space points from the same tower from TkrClusterCol
 *
 * @author Tracy Usher
 *
 * $Header$
 *
*/

#ifndef __TkrTreeTrackFinder_H
#define __TkrTreeTrackFinder_H 1

#include "Event/Recon/TkrRecon/TkrTree.h"
#include "Event/Recon/TkrRecon/TkrBoundBox.h"

#include "TkrVecNodesBuilder.h"

#include "TkrUtil/ITkrGeometrySvc.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "TkrRecon/Track/ITkrFitTool.h"
#include "TkrRecon/Track/IFindTrackHitsTool.h"
#include "src/PatRec/BuildTkrTrack.h"

#include <vector>

// Class definition in the source code file for TkrTreeTrackFinder
class TkrTreePosition;

class TkrTreeTrackFinder 
{
public:
    TkrTreeTrackFinder(IDataProviderSvc*      dataSvc, 
                       ITkrGeometrySvc*       geoSvc,
                       ITkrQueryClustersTool* clusTool,
                       ITkrFitTool*           trackFitTool,
                       IFindTrackHitsTool*    findHitsTool, 
                       Event::TkrClusterCol*  clusterCol);

    ~TkrTreeTrackFinder();

    /// Method to build the tree objects
    int    findTracks(Event::TkrTree* tree, double eventEnergy);

private:

    /// Make the std::set of leaves for a given tree
    int makeLeafSet(Event::TkrVecNodeSet&           leafSet,
                    const Event::TkrNodeSiblingMap* siblingMap);

    /// Used to set the best and next best branch bits after leaf finding
    void setBranchBits(Event::TkrVecNode* node, bool isMainBranch);

    /// Use this to define a look up object for used clusters
    typedef std::set<const Event::TkrCluster*> UsedClusterList;

    /// Recursive routine for building a node sibling map from only the best branch
    void leafToSiblingMap(Event::TkrNodeSiblingMap* siblingMap, 
                          const Event::TkrVecNode*  headNode,
                          UsedClusterList&          usedClusters);

    /// This makes a TkrTrack from the nodes given it in a TkrNodeSiblingMap
    Event::TkrTrack* getTkrTrackFromLeaf(Event::TkrVecNode* leaf, double energy, UsedClusterList& usedClusters);

    /// This makes a TkrTrack, given a start point and direction, by using the kalman hit finder
    Event::TkrTrack* getTkrTrackFromHits(Point  startPoint, Vector startDir, double energy);
    Event::TkrTrack* getTkrTrackFromHits2(const Event::TkrVecPointsLink* firstLink, double energy);

    /// This makes a TkrTrack using a TkrNodeSiblingMap as a guide, depending on nRequiredHits
    Event::TkrTrack* makeTkrTrack(Event::TkrNodeSiblingMap* siblingMap,
                                  Event::TkrFilterParams*   filterParams,
                                  double                    energy        = 1000., 
                                  int                       nRequiredHits = 5);
    Event::TkrTrack* makeTkrTrackFromMean(Event::TkrNodeSiblingMap* siblingMap,
                                          Event::TkrFilterParams*   filterParams,
                                          double                    energy        = 1000., 
                                          int                       nRequiredHits = 5);

    /// Use this to try to set up the tree axis
    typedef std::list<Event::TkrBoundBox*> TkrBoundBoxList;

    void findTreeAxis(Event::TkrNodeSiblingMap* siblingMap, TkrBoundBoxList& bboxList, Point& centroid);

    Event::TkrFilterParams* doMomentsAnalysis(TkrBoundBoxList& bboxList, Point& centroid);

    /// Build the candidate track hit vector which is used to make TkrTracks
    BuildTkrTrack::CandTrackHitVec getCandTrackHitVec(Event::TkrNodeSiblingMap* siblingMap,
                                                      Event::TkrFilterParams*   filterParams);
    BuildTkrTrack::CandTrackHitVec getCandTrackHitVecFromMean(Event::TkrNodeSiblingMap* siblingMap,
                                                              Event::TkrFilterParams*   filterParams);
    BuildTkrTrack::CandTrackHitVec getCandTrackHitVecFromLeaf(Event::TkrVecNode* leaf, UsedClusterList& usedClusters);

    BuildTkrTrack::CandTrackHitVec getCandTrackHitVecFromFirstLink(const Event::TkrVecPointsLink* vecLink);

    /// Use this to handle links that skip layes in the above two methods
    void handleSkippedLayers(const Event::TkrVecPointsLink* vecLink, BuildTkrTrack::CandTrackHitVec& clusVec);

    /// Attempt to not repeat code... 
    void insertVecPointIntoClusterVec(const Event::TkrVecPoint*       vecPoint, 
                                      BuildTkrTrack::CandTrackHitVec& clusVec,
                                      UsedClusterList&                usedClusters);

    /// For calculating the initial position and direction to give to the track
    typedef std::pair<Point, Vector> TkrInitParams;
    TkrInitParams getInitialParams(BuildTkrTrack::CandTrackHitVec& clusVec);

    /// Use this to flag the used clusters as they get used by found tracks
    void flagUsedClusters(UsedClusterList& usedClusters);

    /// Use this to flag all clusters in a given tree (as the last step)
    void flagAllUsedClusters(const Event::TkrTree* tree);

    /// Makes a TkrId
    idents::TkrId makeTkrId(Point& planeHit, int planeId);

    /// Data provider service
    IDataProviderSvc*      m_dataSvc;

    /// Pointer to the local Tracker geometry service
    ITkrGeometrySvc*       m_tkrGeom;

    /// Query Clusters tool
    ITkrQueryClustersTool* m_clusTool;

    /// Track fit tool
    ITkrFitTool*           m_trackFitTool;

    /// Hit finder tool
    IFindTrackHitsTool*    m_findHitsTool;

    /// Pointer to the local cluster collection that we manage
    Event::TkrClusterCol*  m_clusterCol;

    /// Parameter to determine whether to try a composite track
    double                 m_maxChiSqSeg4Composite;

    /// Parameter to control when to allow hit finder to add hits
    double                 m_maxFilterChiSqFctr;

    /// Parameter to control number of shared leading hits
    int                    m_maxSharedLeadingHits;

    /// Parameter used in track hit finding determining maximum gaps between hits
    int                    m_maxGaps;

    /// Parameter used in track hit finding determing maximum consecutive gaps
    int                    m_maxConsecutiveGaps;
};

#endif