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

#include "Event/TopLevel/EventModel.h"
#include "Event/Recon/CalRecon/CalXtalRecData.h"
#include "Event/Recon/CalRecon/CalCluster.h"

#include <CalRecon/ICalReconSvc.h>
#include <CalRecon/ICalClusteringTool.h>
#include "StdClusterInfo.h"
#include "MomentsClusterInfo.h"


// A useful typedef 
typedef  XtalDataList::iterator XtalDataListIterator;

// Define the map between the layer/row index and the crystals
typedef std::map <int, XtalDataList> LyrRow2XtalDataMap;
typedef std::pair<int, XtalDataList> LyrRow2XtalDataPair;


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
    //! Collect CalXtalRecData pointers
    void getXtals();
    //! Distinguish sets of related xtals
    void makeSets(XtalDataListVec& clusters ) ;
    // The driving method for associating crystals in a given layer
    Event::CalXtalRecData* getNearestXtalInSameLayer(XtalDataList& NNvec, Event::CalXtalRecData* xTal);
    // Used by the above driver to associate crystals in the same row
    Event::CalXtalRecData* getNearestXtalInSameRow(XtalDataList& NNvec, Event::CalXtalRecData* xTal);
    // Used by the above driver to find a crystal in the "next" row to use as a seed
    Event::CalXtalRecData* getNearestXtalInNextRow(XtalDataList& NNvec, Event::CalXtalRecData* xTal, int rowInc = 1);
    // Will find the closest crystal in the layer above
    Event::CalXtalRecData* getNearestXtalInDiffLayer(Event::CalXtalRecData* xTal, int layer);
    // uses the Xtal identifier to build a layer/row index
    int  getXtalLayerRow(Event::CalXtalRecData* xTal);
    // remove a crystal from the lists
    void removeXTal(Event::CalXtalRecData* xTal);

    //! Service for basic Cal info
    ICalReconSvc*            m_calReconSvc;

    //! Utility for filling clusters
    ICalClusterFiller*       m_clusterInfo;

    //! Event Service member directly useable by concrete classes.
    IDataProviderSvc*        m_dataSvc;

    // Keep trackk of crystals locally
    XtalDataList             m_xTals;
    LyrRow2XtalDataMap       m_lyrRow2XtalDataMap;
    Event::CalClusterHitTab* m_xTal2ClusTab;

    // Scale factor for searching neighboring rows
    double                   m_nextRowSclFctr;
    double                   m_nextLyrSclFctr;
 } ;


DECLARE_TOOL_FACTORY(CalSimpleClusteringTool) ;

//
// Constructor
//

CalSimpleClusteringTool::CalSimpleClusteringTool(const std::string & type,
                                                 const std::string & name,
                                                 const IInterface * parent )
                                               : AlgTool(type, name, parent),
                                                 m_nextRowSclFctr(0.75),
                                                 m_nextLyrSclFctr(0.5)
{ 
    m_xTals.clear();
    m_lyrRow2XtalDataMap.clear();
    m_xTal2ClusTab = 0;

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

    if(service( "EventDataSvc", m_dataSvc, true ).isFailure()) 
    {
        throw GaudiException("Service [EventDataSvc] not found", name(), sc);
    }

    // Cluster filling utility
    // -- To be replaced with a generic version soon
    m_clusterInfo = new MomentsClusterInfo(m_calReconSvc);

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

    // create a new cluster <--> crystal relational table

    // Create a Xtal to Cluster relations list
    Event::CalClusterHitTabList* xTal2ClusTabList = new Event::CalClusterHitTabList();
    xTal2ClusTabList->clear();

    // Now create the table to use it
    m_xTal2ClusTab = new Event::CalClusterHitTab(xTal2ClusTabList);

    // Register the list in the TDS (which will assume ownership of the list, but not the table)
    if (m_dataSvc->registerObject(EventModel::CalRecon::CalClusterHitTab, xTal2ClusTabList).isFailure())
    {
        throw GaudiException("Unable to register xTal to Cluster table in TDS", name(), StatusCode::FAILURE);
    }

    // "get" the crystals and put in our internal lists
    getXtals();
  
    // find out groups, with at least a zero
    // cluster so that downstream code runs ok
    XtalDataListVec clusters ;
    if (!m_xTals.empty())
    { 
        makeSets(clusters) ; 
    }
    else
    { 
        clusters.push_back(new XtalDataList) ; 
    }

    // Convert the results into CalClusters
    calClusterCol->clear() ;

    for (XtalDataListVec::iterator xTalClusIter = clusters.begin(); xTalClusIter != clusters.end(); xTalClusIter++)
    {
        XtalDataList* xTalClus = *xTalClusIter;

        Event::CalCluster* cluster = m_clusterInfo->fillClusterInfo(xTalClus);
        std::string producerName("CalSimpleClusteringTool/") ;
        producerName += cluster->getProducerName() ;
        cluster->setProducerName(producerName) ;
		cluster->clearStatusBit(Event::CalCluster::ALLXTALS);

        calClusterCol->push_back(cluster);

        // Loop through the xtals to make the relational table (hmmm... this could be done better...)
        for(XtalDataListIterator xTalIter = xTalClus->begin(); xTalIter != xTalClus->end(); xTalIter++)
        {
            Event::CalXtalRecData*   xTal         = *xTalIter;
            Event::CalClusterHitRel* xTal2ClusRel = new Event::CalClusterHitRel(xTal,cluster);

            m_xTal2ClusTab->addRelation(xTal2ClusRel);
        }

        delete xTalClus;
    }

    // Clear the local clusters vector
    clusters.clear();

    // And delete the table (remembering that the list is in the TDS)
    delete m_xTal2ClusTab;
  
    return StatusCode::SUCCESS ;
}

//
// Define a class for the sorting algorithm
// This will order our crystals in our crystal list from highest to lowest energy
//
class CompareXtalEnergy
{
public:
    const bool operator()(const Event::CalXtalRecData* left, const Event::CalXtalRecData* right) const
    {
        if (left->getEnergy() > right->getEnergy()) return true;
        
        return false;
    }
};

//! Collect CalXtalRecData pointers
void CalSimpleClusteringTool::getXtals()
{
    // Make sure we have clean vectors to begin
    m_xTals.clear();
    m_lyrRow2XtalDataMap.clear();

    // Now go through the input collection of xTals
    Event::CalXtalRecCol::const_iterator it ;
    for ( it = m_calReconSvc->getXtalRecs()->begin() ; it != m_calReconSvc->getXtalRecs()->end() ; ++it )
    {
        // get pointer to the reconstructed data for given crystal
	    Event::CalXtalRecData* recData = *it ;

        XtalDataListIterator xTalIter = m_xTals.insert(m_xTals.end(), recData);

        // get a layer/row index to use to associate this crystal
        int layerRow = getXtalLayerRow(recData);

        // Perhaps some paranoia here... but if this is new index then explicitly create the pair in the map
        if (m_lyrRow2XtalDataMap.find(layerRow) == m_lyrRow2XtalDataMap.end())
        {
            m_lyrRow2XtalDataMap.insert(LyrRow2XtalDataPair(layerRow, std::list<Event::CalXtalRecData*>()));
        }

        // Add this crystal to our map
        m_lyrRow2XtalDataMap[layerRow].push_back(recData);
    }

    // Ok, now sort highest to lowest energy
    m_xTals.sort(CompareXtalEnergy());
}

int  CalSimpleClusteringTool::getXtalLayerRow(Event::CalXtalRecData* xTal)
{
    idents::CalXtalId xTalId = xTal->getPackedId();
        
    int curLayer  = xTalId.getLayer();
    int curTower  = xTalId.getTower();
    int row       = xTalId.isX() ? curTower % 4 : curTower / 4;
    int layerRow  = 4 * curLayer + row;

    return layerRow;
}

void CalSimpleClusteringTool::removeXTal(Event::CalXtalRecData* xTal)
{
    // First task is to remove from the map between layer/row and xTals
    // Get the layerRow index
    int layerRow = getXtalLayerRow(xTal);

    // Look this index up in the map to get the vector of xTals
    LyrRow2XtalDataMap::iterator mapIter = m_lyrRow2XtalDataMap.find(layerRow);

    // Make sure we have something (should throw an error if nothing found? 
    if (mapIter != m_lyrRow2XtalDataMap.end())
    {
        // Find the actual xTal in the vector
        std::list<Event::CalXtalRecData*>::iterator xTalDataIter = std::find(mapIter->second.begin(), mapIter->second.end(), xTal);

        // More safety nonsense
        if (xTalDataIter != mapIter->second.end())
        {
            // Actuall remove it
            mapIter->second.erase(xTalDataIter);
        }
    }
    else
    {
        int j = 0;
    }

    // Now remove the crystal from the orginal list
    XtalDataList::iterator xTalDataIter = std::find(m_xTals.begin(), m_xTals.end(), xTal);

    if (xTalDataIter != m_xTals.end())
    {
        m_xTals.erase(xTalDataIter);
    }
    else
    {
        int j = 0;
    }

    return;
}

void CalSimpleClusteringTool::makeSets(XtalDataListVec& clusters)
{
    // The first set is the set of ALL crystals from which we will construct the "uber" cluster
    // Note that this cluster will be added to the end of our list, at the end of this method
    XtalDataList* uber = new XtalDataList(m_xTals);

    // Now make the "proper" clusters
    while (m_xTals.size()>0)
    {
        XtalDataList* cluster = new XtalDataList;

        // The highest energy crystal in our list is the first one. 
        XtalDataList::iterator xTalBestIter = m_xTals.begin();
        Event::CalXtalRecData* xTalBestEne  = *xTalBestIter;

        //Add this crystal to our new cluster list and remove from the old list
        Event::CalXtalRecData* xTal = *xTalBestIter;
        cluster->push_back(xTal);
        removeXTal(xTal);

        //Now build up a list of all xTals in this layer which are neighbors 
        Event::CalXtalRecData* nxtXtal = getNearestXtalInSameLayer(*cluster, xTal);

        //Get the current xTal id
        const idents::CalXtalId xTalId = xTal->getPackedId();

        //Loop "up" (if possible) associating crystals above us
        Event::CalXtalRecData* bestXtal = xTal;
        for(int layer = xTalId.getLayer() - 1; layer >= 0; layer--)
        {
            //Finds the closest xTal in the next layer, if successful removes it from list
            bestXtal = getNearestXtalInDiffLayer(bestXtal, layer);

            //No crystal signifies we are done
            if (bestXtal == 0) break;
            //if (bestXtal == 0) continue;

            //Add to cluster
            removeXTal(bestXtal);
            cluster->push_back(bestXtal);

            //Searches for nearest neigbors, removes from current list when found
            bestXtal = getNearestXtalInSameLayer(*cluster, bestXtal);
        }

        //Loop "down" (if possible) associating xTals below us
        bestXtal = xTal;
        for(int layer = xTalId.getLayer() + 1; layer < m_calReconSvc->getCalNLayers() ; layer++)
        {
            bestXtal = getNearestXtalInDiffLayer(bestXtal, layer);

            if (bestXtal == 0) break;

            removeXTal(bestXtal);
            cluster->push_back(bestXtal);

            bestXtal = getNearestXtalInSameLayer(*cluster, bestXtal);
        }

        clusters.push_back(cluster) ;
    }

    // Now add the uber cluster to the end of our list
    // But only do so if more than one cluster found
    if (clusters.size() > 1) clusters.push_back(uber);
    else                     delete uber;

    // Finished! 
    return;
}

Event::CalXtalRecData* CalSimpleClusteringTool::getNearestXtalInSameLayer(XtalDataList& NNvec, Event::CalXtalRecData* xTal)
{
    Event::CalXtalRecData* bestXtal = getNearestXtalInSameRow(NNvec, xTal);

    // Check in row to "left" in row number
    Event::CalXtalRecData* xTalLower = getNearestXtalInNextRow(NNvec, bestXtal, -1);

    if (xTalLower)
    {
        removeXTal(xTalLower);
        NNvec.push_back(xTalLower);

        xTalLower = getNearestXtalInSameRow(NNvec, xTalLower);
    }

    // Check in row to "right" in row number
    Event::CalXtalRecData* xTalUpper = getNearestXtalInNextRow(NNvec, bestXtal, 1);

    if (xTalUpper)
    {
        removeXTal(xTalUpper);
        NNvec.push_back(xTalUpper);

        xTalUpper = getNearestXtalInSameRow(NNvec, xTalUpper);
    }

    // Now check for "best"
    if (xTalLower && xTalLower->getEnergy() > bestXtal->getEnergy()) bestXtal = xTalLower;
    if (xTalUpper && xTalUpper->getEnergy() > bestXtal->getEnergy()) bestXtal = xTalUpper;

    return bestXtal;
}

Event::CalXtalRecData* CalSimpleClusteringTool::getNearestXtalInSameRow(XtalDataList& NNvec, Event::CalXtalRecData* xTal)
{
    // Return pointer to the best crystal
    Event::CalXtalRecData* bestXtal = xTal;

    // Null input pointer means nothing to do
    if (!bestXtal) return bestXtal;

    //Extract the ID information we need from current crystal
    idents::CalXtalId xTalId = xTal->getPackedId();

    int curTower  = xTalId.getTower();
    int curColumn = xTalId.getColumn();
    int logNum    = xTalId.isX() 
                  ? curColumn + 13 * (curTower / 4) 
                  : curColumn + 13 * (curTower % 4);

    // Retreive the layer/row index for this crystal
    int layerRow = getXtalLayerRow(xTal);

    // Pointer to the "next" found Xtal
    Event::CalXtalRecData* nextXtal = 0;

    // Now find the vector of crystals in this row
    LyrRow2XtalDataMap::iterator lyrRow2XtalIter = m_lyrRow2XtalDataMap.find(layerRow);

    // Make sure we have something (should throw an error if nothing found? 
    if (lyrRow2XtalIter != m_lyrRow2XtalDataMap.end())
    {
        // Get reference to our list
        XtalDataList& xTalDataList = lyrRow2XtalIter->second;
        
        //Loop through input vector of Xtals looking for a nearest neighbor match
        XtalDataListIterator xTalVecIter = xTalDataList.begin();
        while(xTalVecIter != xTalDataList.end())
        {
            Event::CalXtalRecData* xTalTest   = *xTalVecIter;
            idents::CalXtalId      xTalTestId = xTalTest->getPackedId();

            int column    = xTalTestId.getColumn(); 
            int thisTower = xTalTestId.getTower();
            int curLogNum = xTalTestId.isX() 
                          ? column + 13 * (thisTower / 4) 
                          : column + 13 * (thisTower % 4);
            int logDiff   = curLogNum - logNum;

            // Accept xtals that are "near" the current one
            if (abs(logDiff) < 3)
            {
                // Remove from original list
                removeXTal(xTalTest);

                // Add to our neighbor list
                NNvec.push_back(xTalTest);

                // is this the highest energy crystal in this layer?
                if (bestXtal->getEnergy() < xTalTest->getEnergy()) bestXtal = xTalTest;

                // Find this crystals neighbors
                Event::CalXtalRecData* nextXtal = getNearestXtalInSameRow(NNvec, xTalTest);

                // Check if it beats our current best
                if (bestXtal->getEnergy() < nextXtal->getEnergy()) bestXtal = nextXtal;

                // Reset iterator to front beginning of list
                xTalVecIter = xTalDataList.begin();
            }
            else xTalVecIter++;
        }
    }

    //Return pointer to found crystal (if one)
    return bestXtal;
}

Event::CalXtalRecData* CalSimpleClusteringTool::getNearestXtalInNextRow(XtalDataList& NNvec, Event::CalXtalRecData* xTal, int rowInc)
{
    // Return pointer to the best crystal
    Event::CalXtalRecData* bestXtal = 0;

    // Retreive the layer/row index for this crystal
    int layerRow = getXtalLayerRow(xTal);
    int row      = layerRow % 4;

    if (row + rowInc < 0 || row + rowInc > 3) return bestXtal;

    //Extract the ID information we need from current crystal
    idents::CalXtalId xTalId = xTal->getPackedId();

    // And get the log num for this xTal
    int   curTower  = xTalId.getTower();
    int   curColumn = xTalId.getColumn();
    int   logNum    = xTalId.isX() 
                    ? curColumn + 13 * (curTower / 4) 
                    : curColumn + 13 * (curTower % 4);

    // Also need our current location
    Point curPos   = xTal->getPosition();

    // Pointer to the "next" found Xtal
    Event::CalXtalRecData* nextXtal = 0;

    // Now find the vector of crystals in the next row
    LyrRow2XtalDataMap::iterator lyrRow2XtalIter = m_lyrRow2XtalDataMap.find(layerRow + rowInc);

    // Make sure we have something (should throw an error if nothing found? 
    if (lyrRow2XtalIter != m_lyrRow2XtalDataMap.end())
    {
        // Get reference to our list
        XtalDataList& xTalDataList = lyrRow2XtalIter->second;

        //Keep track of the best match
        double  bestDist = m_nextRowSclFctr * m_calReconSvc->getCalCsILength();  
        
        //Loop through input vector of Xtals looking for a nearest neighbor match
        XtalDataListIterator xTalVecIter = xTalDataList.begin();
        while(xTalVecIter != xTalDataList.end())
        {
            Event::CalXtalRecData* xTalTest   = *xTalVecIter++;
            idents::CalXtalId      xTalTestId = xTalTest->getPackedId();

            int column    = xTalTestId.getColumn(); 
            int thisTower = xTalTestId.getTower();
            int curLogNum = xTalTestId.isX() 
                          ? column + 13 * (thisTower / 4) 
                          : column + 13 * (thisTower % 4);
            int logDiff   = curLogNum - logNum;

            //Compute distance to this xTal
            Vector distVec = curPos - xTalTest->getPosition();
            double dist = distVec.magnitude() / m_calReconSvc->getCalCsIHeight() ;

            // Accept xtals that are "near" the current one
            if (abs(logDiff) < 3 && dist < bestDist)
            {
                bestXtal = xTalTest;
            }
        }
    }

    //Return pointer to found crystal (if one)
    return bestXtal;
}

Event::CalXtalRecData* CalSimpleClusteringTool::getNearestXtalInDiffLayer(Event::CalXtalRecData* xTal, int layer)
{
    Event::CalXtalRecData* newXtal = 0;

    //Extract the ID information we need from current crystal, also get position
    idents::CalXtalId xTalId   = xTal->getPackedId();
    int               curTower = xTalId.getTower();
    Point             curPos   = xTal->getPosition();

    // Build a layer/row index for the desired layer assuming the tower of the xTal passed in
    int layerRow = layer % 2 == 0 ? 4 * layer + (curTower % 4) : 4 * layer + (curTower / 4);

    // Use this to look up the crystals in the next layer
    LyrRow2XtalDataMap::iterator lyrRow2XtalIter = m_lyrRow2XtalDataMap.find(layerRow);

    // Make sure we have something (should throw an error if nothing found? 
    if (lyrRow2XtalIter != m_lyrRow2XtalDataMap.end())
    {
        //Keep track of the best match
        double  bestDist = m_nextLyrSclFctr * m_calReconSvc->getCalCsILength();   

        //Loop through input vector of Xtals looking for a nearest neighbor match
        std::list<Event::CalXtalRecData*>::iterator xTalVecIter = lyrRow2XtalIter->second.begin();
        while(xTalVecIter != lyrRow2XtalIter->second.end())
        {
            Event::CalXtalRecData*  xTalTest   = *xTalVecIter++;
            const idents::CalXtalId xTalTestId = xTalTest->getPackedId();

            //Compute distance to this xTal
            Vector distVec = curPos - xTalTest->getPosition();
            double dist = distVec.magnitude() / m_calReconSvc->getCalCsIHeight() ;

            if (dist < bestDist)
            {
                bestDist = dist;
                newXtal  = xTalTest;
            }
        }
    }

    //Return out find...
    return newXtal;
}

