/// @file TkrPoints.cxx
/**
* @brief Provides X-Y space points from the same tower from TkrClusterCol
*
* 1-Dec-2001
*  This class provides an easy way to cycle over allowed XY pairs
*  in a given GLAST paired layer.  The ordering can be either 
*  combinatoric or based on nearest, next nearest, etc. to a given
*  point in that plane

* 11/01/2004
*  Rewritten for new scheme
*  LSR
*
* @authors Bill Atwood, Brian Algood
*
* $Header$
*
*/

#include "src/Utilities/TkrPoints.h"
#include <iostream>

#include <algorithm>

TkrPoints::TkrPoints(int layer, ITkrQueryClustersTool* clusTool)
{
    m_layer    = layer;
    m_clusTool = clusTool;
    m_refPoint = Point(0, 0, 0);
    m_maxDist2 = _big*_big;
    m_sorting  = false;
    ini();
}


TkrPoints::TkrPoints(int layer, ITkrQueryClustersTool* clusTool, 
                     Point refPoint, double maxDistance)
{
    m_layer    = layer;
    m_clusTool = clusTool;
    m_refPoint = refPoint;
    m_maxDist2 = maxDistance*maxDistance;
    m_sorting  = true;
    ini();
    // sort by closest to refPoint (see header file for predicate)
    std::sort(this->begin(), this->end(), sortByClosest(m_refPoint));    
}

void TkrPoints::ini()
{
    this->clear();

    Event::TkrClusterVec xHitList 
        = m_clusTool->getClusters(idents::TkrId::eMeasureX,m_layer);

    Event::TkrClusterVec yHitList 
        = m_clusTool->getClusters(idents::TkrId::eMeasureY,m_layer);

    Event::TkrClusterVecConItr itX = xHitList.begin();
    for (; itX!=xHitList.end(); ++itX) {
        const Event::TkrCluster* clX = *itX;
        Event::TkrClusterVecConItr itY = yHitList.begin();
        for (; itY!=yHitList.end(); ++itY) {
            const Event::TkrCluster* clY = *itY;
            if(clX->tower()!=clY->tower()) continue;
            TkrPoint* pt = new TkrPoint(m_layer, clX, clY);
            // in sorting mode, throw out far-away points (defaults to keep all)
            if(m_sorting) {
                if(pt->getDistanceSquaredTo(m_refPoint)>m_maxDist2) {
                    delete pt;
                    continue;
                }
            }
            this->push_back(pt);
        }
    }
}
   
TkrPoint* TkrPoints::getNearestPointOutside(Point x0, double & dist_min) const
{
    // Searches out the nearest space point to x0 which lies
    // outside a distance d
    // returns distance to this point, negative if no point is found

    double dist_best2 = _big*_big;
    double dist_min2 = dist_min*dist_min;
    bool found = false;
    TkrPoint* ret = 0;

    TkrPointListConItr iPt = this->begin();
    for (; iPt!=this->end(); ++iPt) {
        TkrPoint* pt = (*iPt);
        if (pt->flagged()) continue;
        double dist2 = pt->getDistanceSquaredTo(x0);
        if(dist2 < dist_min2 || dist2 > dist_best2) continue; 
        found = true;
        dist_best2 = dist2;
        ret = pt;
    }

    dist_min = (found ? sqrt(dist_best2) : -1.0 );
    return ret;
}
