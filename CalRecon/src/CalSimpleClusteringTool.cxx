/**
 * @class CalSimpleClusteringTool
 *
 * @brief Implements a Gaudi Tool for performing very simple clustering in the Cal 
 *
 * @author Tracy Usher
 *
 * File and Version Information:
 *      $Header$
 */

// Tool and Gaudi related stuff
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/GaudiException.h" 
#include "GaudiKernel/IParticlePropertySvc.h"
#include "GaudiKernel/DeclareFactoryEntries.h"
#include "CalSimpleClusteringTool.h"

DECLARE_TOOL_FACTORY(CalSimpleClusteringTool) ;

//
// Feeds Combo pattern recognition tracks to Kalman Filter
//

CalSimpleClusteringTool::CalSimpleClusteringTool
 ( const std::string & type,
   const std::string & name,
   const IInterface * parent )
 : CalClusteringTool(type, name, parent)
 { declareInterface<CalIClusteringTool>(this) ; }

// 
// Cleanup memory on exit
//
CalSimpleClusteringTool::~CalSimpleClusteringTool()
 {}
 
/// This finds the next highest energy cluster in a vector of CalXtalRecData pointers
CalClusteringTool::xTalDataVec CalSimpleClusteringTool::nextXtalsSet(xTalDataVec& xTalVec)
{
    xTalDataVec cluster;
    cluster.clear();

    //Start by finding the highest energy crystal in the vector
    double bestE = 0.;
    xTalDataVec::iterator bestIter;
    xTalDataVec::iterator xTalVecIter;
    for(xTalVecIter = xTalVec.begin(); xTalVecIter != xTalVec.end(); xTalVecIter++)
    {
        Event::CalXtalRecData* xTalRec = *xTalVecIter;

        if (xTalRec->getEnergy() > bestE)
        {
            bestIter = xTalVecIter;
            bestE    = xTalRec->getEnergy();
        }
    }

    //Add this crystal to our new cluster list and remove from the old list
    Event::CalXtalRecData* xTal = *bestIter;
    cluster.push_back(*bestIter);
    xTalVec.erase(bestIter);

    //Now build up a list of all xTals in this layer which are neighbors 
    Event::CalXtalRecData* nxtXtal = getNearestXtalInSameLayer(xTalVec, cluster, xTal);

    if (nxtXtal) cluster.push_back(nxtXtal);

    //Do again to pick up neighbor on the other side (we need a better algorithm here...) 
    nxtXtal = getNearestXtalInSameLayer(xTalVec, cluster, xTal);

    if (nxtXtal) cluster.push_back(nxtXtal);

    //Get the current xTal id
    const idents::CalXtalId xTalId = xTal->getPackedId();

    //Loop "up" (if possible) associating crystals above us
    Event::CalXtalRecData* bestXtal = xTal;
    for(int layer = xTalId.getLayer() - 1; layer >= 0; layer--)
    {
        //Finds the closest xTal in the next layer, if successful removes it from list
        bestXtal = getNearestXtalInDiffLayer(xTalVec, bestXtal, layer);

        //No crystal signifies we are done
        if (bestXtal == 0) break;
        //if (bestXtal == 0) continue;

        //Add to cluster
        cluster.push_back(bestXtal);

        //Searches for nearest neigbors, removes from current list when found
        nxtXtal = getNearestXtalInSameLayer(xTalVec, cluster, bestXtal);
        if (nxtXtal) cluster.push_back(nxtXtal);
        nxtXtal = getNearestXtalInSameLayer(xTalVec, cluster, bestXtal);
        if (nxtXtal) cluster.push_back(nxtXtal);
    }

    //Loop "down" (if possible) associating xTals below us
    bestXtal = xTal;
    for(int layer = xTalId.getLayer() + 1; layer < m_CalnLayers; layer++)
    {
        bestXtal = getNearestXtalInDiffLayer(xTalVec, bestXtal, layer);

        if (bestXtal == 0) break;
        //if (bestXtal == 0) continue;

        cluster.push_back(bestXtal);

        nxtXtal = getNearestXtalInSameLayer(xTalVec, cluster, bestXtal);
        if (nxtXtal) cluster.push_back(nxtXtal);
        nxtXtal = getNearestXtalInSameLayer(xTalVec, cluster, bestXtal);
        if (nxtXtal) cluster.push_back(nxtXtal);
    }

    // Finished! 
    return cluster;
}

CalClusteringTool::xTalDataVec CalSimpleClusteringTool::getXtalsInLayer(xTalDataVec& xTalVec, Event::CalXtalRecData* xTal)
{
    xTalDataVec newVec;
    newVec.clear();

    //Extract the ID information we need from current crystal
    const idents::CalXtalId xTalId   = xTal->getPackedId();
    int                     curTower = xTalId.getTower();
    int                     curLayer = xTalId.getLayer();

    //Loop through input vector of Xtals looking for a match to this layer
    xTalDataVec::iterator xTalVecIter;
    for(xTalVecIter = xTalVec.begin(); xTalVecIter != xTalVec.end(); xTalVecIter++)
    {
        Event::CalXtalRecData* xTalRec = *xTalVecIter;

        const idents::CalXtalId xTalId  = xTalRec->getPackedId();

        //If same tower and layer then store away
        if (xTalId.getTower() == curTower && xTalId.getLayer() == curLayer) newVec.push_back(xTalRec);
    }

    return newVec;
}

Event::CalXtalRecData* CalSimpleClusteringTool::getNearestXtalInSameLayer(xTalDataVec& xTalVec, xTalDataVec& NNvec, Event::CalXtalRecData* xTal)
{
    //Extract the ID information we need from current crystal
    const idents::CalXtalId xTalId    = xTal->getPackedId();
    int                     curTower  = xTalId.getTower();
    int                     curLayer  = xTalId.getLayer();
    int                     curColumn = xTalId.getColumn();

    //Loop through input vector of Xtals looking for a nearest neighbor match
    xTalDataVec::iterator xTalVecIter;
    for(xTalVecIter = xTalVec.begin(); xTalVecIter != xTalVec.end(); xTalVecIter++)
    {
        Event::CalXtalRecData* xTalRec  = *xTalVecIter;
        const idents::CalXtalId xTalId  = xTalRec->getPackedId();
        int                     column  = xTalId.getColumn(); 
        int                     colDiff = curColumn - column;

        //Only accept Xtal if it is exactly next to the current one
        if (curTower == xTalId.getTower() && curLayer == xTalId.getLayer() && abs(colDiff) < 2)
        {
            //Remove this Xtal from the current list
            xTalVec.erase(xTalVecIter);

            //Look for the nearest neighbor to this Xtal
            Event::CalXtalRecData* nextXtal = getNearestXtalInSameLayer(xTalVec, NNvec, xTalRec);

            //If one found, add to the Nearest Neighbor list
            if (nextXtal) NNvec.push_back(nextXtal);

            //return the current Xtal
            return xTalRec;
        }
    }

    //If we got here then nothing found, return a null pointer
    return 0;
}

Event::CalXtalRecData* CalSimpleClusteringTool::getNearestXtalInDiffLayer(xTalDataVec& xTalVec, Event::CalXtalRecData* xTal, int layer)
{
    Event::CalXtalRecData* newXtal = 0;

    //Extract the ID information we need from current crystal, also get position
    const idents::CalXtalId xTalId   = xTal->getPackedId();
    int                     curTower = xTalId.getTower();
    Point                   curPos   = xTal->getPosition();
    //bool                    isXlyr   = layer % 2 == 0;

    //Keep track of the best match
    xTalDataVec::iterator bestXtalIter = xTalVec.end();
    double                bestDist     = 10.;

    //Loop through input vector of Xtals looking for a nearest neighbor match
    xTalDataVec::iterator xTalVecIter;
    for(xTalVecIter = xTalVec.begin(); xTalVecIter != xTalVec.end(); xTalVecIter++)
    {
        Event::CalXtalRecData*  xTalRec = *xTalVecIter;
        const idents::CalXtalId xTalId  = xTalRec->getPackedId();

        //Only accept Xtal if it is exactly next to the current one
        if (curTower == xTalId.getTower() && layer == xTalId.getLayer())
        {
            //Compute distance to this xTal
            Vector distVec = curPos - xTalRec->getPosition();
            double dist    = distVec.magnitude() / m_CsIHeight;

            if (dist < bestDist)
            {
                bestDist     = dist;
                bestXtalIter = xTalVecIter;
            }
        }
    }

    //Did we get an acceptable crystal?
    if (bestXtalIter != xTalVec.end())
    {
        newXtal = *bestXtalIter;
        xTalVec.erase(bestXtalIter);
    }

    //Return out find...
    return newXtal;
}

