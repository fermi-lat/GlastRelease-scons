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
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/GaudiException.h" 
#include "GaudiKernel/DeclareFactoryEntries.h"

#include "Event/Recon/CalRecon/CalXtalRecData.h"
#include "Event/Recon/CalRecon/CalCluster.h"

#include "Utilities/ICalReconSvc.h"
#include "ICalClusteringTool.h"
#include "StdClusterInfo.h"

class CalSimpleClusteringTool : public AlgTool, virtual public ICalClusteringTool
{
public :
  
    /// Standard Gaudi Tool interface constructor
    CalSimpleClusteringTool(const std::string& type,
                            const std::string& name,
                            const IInterface* parent );

    virtual ~CalSimpleClusteringTool() {};
    
	/// @brief Intialization of the tool
    virtual StatusCode initialize() ;

    /// @brief Default cluster finding framework
    virtual StatusCode findClusters(Event::CalClusterCol* calClusterCol) ;

private:
    //! Service for basic Cal info
    ICalReconSvc*     m_calReconSvc;

    //! Utility for filling clusters
    IClusterFiller*   m_clusterInfo;

    //! Collect CalXtalRecData pointers
    void getXtals(XtalDataVec& xtals);

    //! Distinguish sets of related xtals
    void makeSets( const XtalDataVec& xtals, XtalDataVecVec& clusters ) ;
  
    XtalDataVec            getXtalsInLayer(XtalDataVec& xTalVec, Event::CalXtalRecData* xTal);
    Event::CalXtalRecData* getNearestXtalInSameLayer(XtalDataVec& xTalVec, XtalDataVec& NNvec, Event::CalXtalRecData* xTal);
    Event::CalXtalRecData* getNearestXtalInDiffLayer(XtalDataVec& xTalVec, Event::CalXtalRecData* xTal, int layer);

 } ;


DECLARE_TOOL_FACTORY(CalSimpleClusteringTool) ;

//
// Constructor
//

CalSimpleClusteringTool::CalSimpleClusteringTool(const std::string & type,
                                                 const std::string & name,
                                                 const IInterface * parent )
                                               : AlgTool(type, name, parent)
{ 
    declareInterface<ICalClusteringTool>(this) ; 
}
    
StatusCode CalSimpleClusteringTool::initialize()
{
    MsgStream log(msgSvc(),name()) ;
    StatusCode sc = StatusCode::SUCCESS;

    if ((sc = service("CalReconSvc",m_calReconSvc,true)).isFailure())
    {
        throw GaudiException("Service [CalReconSvc] not found", name(), sc);
    }

    // Cluster filling utility
    // -- To be replaced with a generic version soon
    m_clusterInfo = new StdClusterInfo(m_calReconSvc);

    return StatusCode::SUCCESS ;
}

StatusCode CalSimpleClusteringTool::findClusters(Event::CalClusterCol* calClusterCol)
{
    //Purpose and method:
    //
    //   This function performs the calorimeter cluster reconstruction.
    //   The main actions are:
    //      - calculate energy sum
    //                  energy per layer
    //                  average position per layer
    //                  quadratic spread per layer
    //      - fit the particle direction using Fit_Direction() function
    //      - store all calculated quantities in CalCluster objects
    // 
    // TDS input: CalXtalRecCol
    // TDS output: CalClustersCol
    // prepare the initital set of xtals
    XtalDataVec xtals ;
    getXtals(xtals) ;
  
    // find out groups, with at least a zero
    // cluster so that downstream code runs ok
    XtalDataVecVec clusters ;
    if (xtals.size()>0)
    { 
        makeSets(xtals,clusters) ; 
    }
    else
    { 
        clusters.push_back(new XtalDataVec) ; 
    }

    // Convert the results into CalClusters
    for (XtalDataVecVec::iterator xTalClusIter = clusters.begin(); xTalClusIter != clusters.end(); xTalClusIter++)
    {
        XtalDataVec* xTalClus = *xTalClusIter;

        Event::CalCluster* cluster = m_clusterInfo->fillClusterInfo(xTalClus);

        calClusterCol->push_back(cluster);

        delete xTalClus;
    }

    clusters.clear();
  
    return StatusCode::SUCCESS ;
}

//! Collect CalXtalRecData pointers
void CalSimpleClusteringTool::getXtals(XtalDataVec& xtals)
{
    xtals.clear();
    Event::CalXtalRecCol::const_iterator it ;
    for ( it = m_calReconSvc->getXtalRecs()->begin() ; it != m_calReconSvc->getXtalRecs()->end() ; ++it )
    {
        // get pointer to the reconstructed data for given crystal
	    Event::CalXtalRecData * recData = *it ;
        xtals.push_back(recData) ;
    }
}

void CalSimpleClusteringTool::makeSets(const XtalDataVec& xtals, XtalDataVecVec& clusters)
{
    XtalDataVec xTalVec(xtals) ;
    while (xTalVec.size()>0)
    {
        XtalDataVec * cluster = new XtalDataVec ;

        //Start by finding the highest energy crystal in the vector
        double bestE = 0.;
        XtalDataVec::iterator bestIter;
        XtalDataVec::iterator xTalVecIter;
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
        cluster->push_back(*bestIter);
        xTalVec.erase(bestIter);

        //Now build up a list of all xTals in this layer which are neighbors 
        Event::CalXtalRecData* nxtXtal = getNearestXtalInSameLayer(xTalVec, *cluster, xTal);

        if (nxtXtal) cluster->push_back(nxtXtal);

        //Do again to pick up neighbor on the other side (we need a better algorithm here...) 
        nxtXtal = getNearestXtalInSameLayer(xTalVec, *cluster, xTal);

        if (nxtXtal) cluster->push_back(nxtXtal);

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
            cluster->push_back(bestXtal);

            //Searches for nearest neigbors, removes from current list when found
            nxtXtal = getNearestXtalInSameLayer(xTalVec, *cluster, bestXtal);
            if (nxtXtal) cluster->push_back(nxtXtal);
            nxtXtal = getNearestXtalInSameLayer(xTalVec, *cluster, bestXtal);
            if (nxtXtal) cluster->push_back(nxtXtal);
        }

        //Loop "down" (if possible) associating xTals below us
        bestXtal = xTal;
        for(int layer = xTalId.getLayer() + 1; layer < m_calReconSvc->getCalNLayers() ; layer++)
        {
            bestXtal = getNearestXtalInDiffLayer(xTalVec, bestXtal, layer);

            if (bestXtal == 0) break;
            //if (bestXtal == 0) continue;

            cluster->push_back(bestXtal);

            nxtXtal = getNearestXtalInSameLayer(xTalVec, *cluster, bestXtal);
            if (nxtXtal) cluster->push_back(nxtXtal);
            nxtXtal = getNearestXtalInSameLayer(xTalVec, *cluster, bestXtal);
            if (nxtXtal) cluster->push_back(nxtXtal);
        }

        clusters.push_back(cluster) ;
    }
    // Finished! 

    return;
}

XtalDataVec CalSimpleClusteringTool::getXtalsInLayer(XtalDataVec& xTalVec, Event::CalXtalRecData* xTal)
{
    XtalDataVec newVec;
    newVec.clear();

    //Extract the ID information we need from current crystal
    const idents::CalXtalId xTalId   = xTal->getPackedId();
    int                     curTower = xTalId.getTower();
    int                     curLayer = xTalId.getLayer();

    //Loop through input vector of Xtals looking for a match to this layer
    XtalDataVec::iterator xTalVecIter;
    for(xTalVecIter = xTalVec.begin(); xTalVecIter != xTalVec.end(); xTalVecIter++)
    {
        Event::CalXtalRecData* xTalRec = *xTalVecIter;

        const idents::CalXtalId xTalId  = xTalRec->getPackedId();

        //If same tower and layer then store away
        if (xTalId.getTower() == curTower && xTalId.getLayer() == curLayer) newVec.push_back(xTalRec);
    }

    return newVec;
}

Event::CalXtalRecData* CalSimpleClusteringTool::getNearestXtalInSameLayer(XtalDataVec& xTalVec, XtalDataVec& NNvec, Event::CalXtalRecData* xTal)
{
    //Extract the ID information we need from current crystal
    const idents::CalXtalId xTalId    = xTal->getPackedId();
    int                     curTower  = xTalId.getTower();
    int                     curLayer  = xTalId.getLayer();
    int                     curColumn = xTalId.getColumn();

    //Loop through input vector of Xtals looking for a nearest neighbor match
    XtalDataVec::iterator xTalVecIter;
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

Event::CalXtalRecData* CalSimpleClusteringTool::getNearestXtalInDiffLayer(XtalDataVec& xTalVec, Event::CalXtalRecData* xTal, int layer)
{
    Event::CalXtalRecData* newXtal = 0;

    //Extract the ID information we need from current crystal, also get position
    const idents::CalXtalId xTalId   = xTal->getPackedId();
    int                     curTower = xTalId.getTower();
    Point                   curPos   = xTal->getPosition();
    //bool                    isXlyr   = layer % 2 == 0;

    //Keep track of the best match
    XtalDataVec::iterator bestXtalIter = xTalVec.end();
    double                bestDist     = 10.;

    //Loop through input vector of Xtals looking for a nearest neighbor match
    XtalDataVec::iterator xTalVecIter;
    for(xTalVecIter = xTalVec.begin(); xTalVecIter != xTalVec.end(); xTalVecIter++)
    {
        Event::CalXtalRecData*  xTalRec = *xTalVecIter;
        const idents::CalXtalId xTalId  = xTalRec->getPackedId();

        //Only accept Xtal if it is exactly next to the current one
        if (curTower == xTalId.getTower() && layer == xTalId.getLayer())
        {
            //Compute distance to this xTal
            Vector distVec = curPos - xTalRec->getPosition();
            double dist = distVec.magnitude() / m_calReconSvc->getCalCsIHeight() ;

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

