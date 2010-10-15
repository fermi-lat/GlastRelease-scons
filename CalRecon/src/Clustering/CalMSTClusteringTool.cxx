/**
 * @class CalMSTClusteringTool
 *
 * @brief Implements a Gaudi Tool for performing MST based clustering in the Cal 
 *
 * @author carmelo
 *
 */

// Tool and Gaudi related stuff
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/GaudiException.h" 
#include "GaudiKernel/DeclareFactoryEntries.h"

#include "Event/TopLevel/EventModel.h"
#include "Event/Recon/CalRecon/CalXtalRecData.h"
#include "Event/Recon/CalRecon/CalCluster.h"
#include "Event/Recon/CalRecon/CalMSTreeParams.h"

#include <CalRecon/ICalReconSvc.h>
#include <CalRecon/ICalClusteringTool.h>
#include "StdClusterInfo.h"
#include "MomentsClusterInfo.h"

#include "TMath.h"

// A useful typedef 
typedef  XtalDataList::iterator XtalDataListIterator;


//
// Define a class for MSTEdge
// An Edge is a weight and two xtals (connected by the line with a weight)
//
class MSTEdge
{

public:
  MSTEdge( )
    :node1(0),node2(0),weight(0.){;};
  MSTEdge(Event::CalXtalRecData& cry1, Event::CalXtalRecData& cry2, double w)
    :node1(&cry1),node2(&cry2),weight(w){;}
  MSTEdge(const MSTEdge& other)
    :node1(other.node1),node2(other.node2),weight(other.weight){;}
  ~MSTEdge( )  { };
  
  void setEdge(Event::CalXtalRecData*, Event::CalXtalRecData*, double);
  Event::CalXtalRecData* getNode1(){return node1;};
  Event::CalXtalRecData* getNode2(){return node2;};
  double getWeight(){return weight;};
  // stupid-useless function that I use for debug
  void printSumEnergy() {std::cout << " Ene1 = " << node1->getEnergy() << " Ene2 = " << node2->getEnergy() << std::endl;};

private:
  Event::CalXtalRecData* node1;
  Event::CalXtalRecData* node2;
  double weight;

};

void MSTEdge::setEdge(Event::CalXtalRecData* _node1, Event::CalXtalRecData* _node2, double _weight)
{
  node1 = _node1;
  node2 = _node2;
  weight = _weight;
}
//
// Define a class for MST Tree: 
// a Tree is a list of Edges + something that may be useful later
// a Tree is a cluster, i.e. a cluster is made from the xtals contained in a tree
//
class MSTTree
{

public:
  MSTTree( ) ;
  ~MSTTree( )  {};
  
  void addEdge(MSTEdge &);
  void addNode(Event::CalXtalRecData* _node) ; // implemented using a map
  int size () {return (m_edges.size());}
  void clear() {m_edges.clear(); nodemap.clear();}; // useful reset
  std::list<MSTEdge*> getEdges() {return m_edges;}; 
  std::map<int, Event::CalXtalRecData*> getNodeMap() {return nodemap;};

  ///--------------------------------------------------- 
  ///--- Get methods
  /// Retrieve the total enegy
  double getTotalEnergy() {return m_totalEnergy;}
  /// Retrieve the energy in the crystal with the maximum enegy
  double getMaxEnergy()   {return m_maxXtalEnergy;}  
  /// Retrieve the number of edges
  int getNumEdges()                  {return (m_edges.size());}
  /// Retrieve the minimum edge length
  double  getMinEdgeLength()	     {return m_minEdgeLength;}
  /// Retrieve the minimum edge length
  double  getMaxEdgeLength()	     {return m_maxEdgeLength;}
  /// Retrieve the maximum edge length
  double  getMeanEdgeLength()	     {return m_meanEdgeLength;}
  /// Retrieve the average edge length
  double  getMeanEdgeLengthTrunc()    {return m_meanEdgeLengthTrunc;}
  /// Retrieve the RMS of edges length
  double  getRmsEdgeLength()	     {return m_rmsEdgeLength;}
  /// Retrieve the RMS of edges length  after truncation
  double  getRmsEdgeLengthTrunc()    {return m_rmsEdgeLengthTrunc;}
  /// Check if stats are available
  double getStatsBit() {return m_statsAvailable;}

  // Other methods
  void printNodes();
  void evalStats();
  
  
private:
  // Internal methods
  std::list<MSTEdge*> m_edges;
  std::map<int, Event::CalXtalRecData*> nodemap;
  int  getXtalUniqueId(Event::CalXtalRecData *);
  bool m_statsAvailable;
  // Members for MSTree properties -- length means weight in term of MSTree
  // Note that numberOfEdges is directly available via m_edges.size()
  double m_totalEnergy;
  double m_maxXtalEnergy;
  
  double m_minEdgeLength;
  double m_maxEdgeLength;
  double m_meanEdgeLength;
  double m_meanEdgeLengthTrunc;
  double m_rmsEdgeLength;
  double m_rmsEdgeLengthTrunc;
  
  
};

MSTTree::MSTTree () 
{
  
  m_edges.clear(); // make sure with know the initial state.
  m_totalEnergy   = 0.;
  m_maxXtalEnergy = 0.;
  m_minEdgeLength   = 0;
  m_maxEdgeLength   = 0;
  m_meanEdgeLength  = 0;
  m_meanEdgeLengthTrunc=0;
  m_rmsEdgeLength   = 0;
  m_rmsEdgeLengthTrunc =0;
  m_statsAvailable=false;
  nodemap.clear();
}

void MSTTree::addEdge(MSTEdge &_edg) 
{
  m_edges.push_back(&_edg);
  // add nodes
  addNode(m_edges.back()->getNode1()); 
  addNode(m_edges.back()->getNode2()); 
}

void MSTTree::addNode(Event::CalXtalRecData *_node) 
{
  nodemap[getXtalUniqueId(_node)] = _node;
}

// collect some statistics for sorting, print info etc...
void MSTTree::evalStats()
{
  // characterize first the energy, loop over the cristals
  m_totalEnergy = 0.;
  m_maxXtalEnergy = 0.;  
  std::map<int, Event::CalXtalRecData*>::iterator it;
  for ( it=nodemap.begin() ; it != nodemap.end(); it++ )
    {
      double xtalEnergy = (*it).second->getEnergy();
      m_totalEnergy += xtalEnergy;
      if (xtalEnergy>=m_maxXtalEnergy) m_maxXtalEnergy = xtalEnergy;
    }

  // now get the property of the edges loop over the edges
  // but do so only if there is at least one edge !
  // if there is no edge, i.e. it's a tree/cluster with one cristal
  // then from the constructor the parameters should be already set to 0
  if(getNumEdges()>0)
    { 
    std::list<MSTEdge*>::iterator itedge;
    double length=0.;
    for ( itedge=m_edges.begin() ; itedge!=m_edges.end(); itedge++ )
      {
    	length=(*itedge)->getWeight();
    	if (length>=m_maxEdgeLength) m_maxEdgeLength = length;
    	if (length<=m_minEdgeLength) m_minEdgeLength = length;
    	m_meanEdgeLength+=length;
    	m_rmsEdgeLength+=length*length;
      }
    m_meanEdgeLength/=getNumEdges();
    m_rmsEdgeLength/=getNumEdges();
    m_rmsEdgeLength=sqrt(m_rmsEdgeLength);

    // now get properties truncating outliers at 3*RMS
    int trcnt=0;
    double maxTrunc=m_meanEdgeLength+3*m_rmsEdgeLength;
    for ( itedge=m_edges.begin() ; itedge!=m_edges.end(); itedge++ )
      {
    	if(length<=maxTrunc)
    	  {
    	  trcnt+=1;
    	  m_meanEdgeLengthTrunc+=length;
    	  m_rmsEdgeLengthTrunc+=length*length;    
    	  }
      }
    m_meanEdgeLengthTrunc/=trcnt; 
    m_rmsEdgeLengthTrunc/=trcnt; 
    m_rmsEdgeLengthTrunc=sqrt(m_rmsEdgeLengthTrunc);
    }
  
  // Set statistics bit
  m_statsAvailable=true;
}

void MSTTree::printNodes() 
{
  // show tree content:
  std::map<int,  Event::CalXtalRecData*>::iterator it;
  for ( it=nodemap.begin() ; it != nodemap.end(); it++ )
    std::cout << "--------> MAP: "<< (*it).first << " => " << (*it).second->getEnergy() << " MeV"<< std::endl;

}

// a unique Id to avoid duplication in the nodemap
int  MSTTree::getXtalUniqueId(Event::CalXtalRecData * xTal)
{
  idents::CalXtalId xTalId = xTal->getPackedId();
  // is this packedId unique?
  return xTalId.getPackedId();
}


// MST tree comparison based on total energy
bool compare_total_energy (MSTTree first, MSTTree second)
{
  
  if (first.getTotalEnergy()>second.getTotalEnergy()) return true;
  else return false;
}

// end of my clustering classes

class CalMSTClusteringTool : public AlgTool, virtual public ICalClusteringTool
{
public :
  
    /// Standard Gaudi Tool interface constructor
    CalMSTClusteringTool(const std::string& type,
			 const std::string& name,
			 const IInterface* parent );

    virtual ~CalMSTClusteringTool() {};
    
    /// @brief Intialization of the tool
    virtual StatusCode initialize() ;

    /// @brief Default cluster finding framework
    virtual StatusCode findClusters(Event::CalClusterCol* calClusterCol) ;

private:
  
  // calculate the weight between two xtals
  double xtalsWeight(Event::CalXtalRecData* xTal1, Event::CalXtalRecData* xTal2 );
  // uses the Xtal identifier to build a layer/row index
  //int  getXtalLayerRow(Event::CalXtalRecData* xTal);
  // return a energy-dependent max weigth
  double getWeightThreshold(double energy);
  
  //! Service for basic Cal info
  ICalReconSvc*      m_calReconSvc;
  
  //! Event Service member directly useable by concrete classes.
  IDataProviderSvc*  m_dataSvc;
  
  //! Utility for filling clusters
  ICalClusterFiller* m_clusterInfo;

  // Keep trackk of crystals locally
  XtalDataList             m_xTals;
  XtalDataList             m_xTals_setA;
  Event::CalClusterHitTab* m_xTal2ClusTab;

  // MST Ubertree and list of all edges
  MSTTree m_uberTree;
  std::list<MSTEdge>  m_uberEdges; // CS: maybe I don't need this one...

  // Clusters are based on this list of trees
  std::list<MSTTree> m_clusterTree;

  // single value overwritten by model
  float    m_maxEdgeWeight;
  // above that we do not grow a tree
  int      m_maxNumXtals;
  // Parameters to separate trees
  float    m_maxEdgeWeightModel_thrLE;     // (mm) @ 10 MeV
  float    m_maxEdgeWeightModel_thrPivEne; // (MeV)
  float    m_maxEdgeWeightModel_thrHE;     // (mm) above the pivot

} ;

DECLARE_TOOL_FACTORY(CalMSTClusteringTool) ;

CalMSTClusteringTool::CalMSTClusteringTool(const std::string & type, 
                                                 const std::string & name, 
                                                 const IInterface* parent)
  : AlgTool(type,name,parent),
    m_maxEdgeWeight(300.) // Overwritten by jobOption
{ 

  m_xTals.clear();
  m_xTals_setA.clear();
  m_clusterTree.clear();
  m_uberEdges.clear();
  declareInterface<ICalClusteringTool>(this) ; 
  
  // jobOptions declaration
  declareProperty ("m_maxNumXtals"                 , m_maxEdgeWeight               = 1536.);
  declareProperty ("m_maxEdgeWeightModel_thrLE"    , m_maxEdgeWeightModel_thrLE    = 500. );
  declareProperty ("m_maxEdgeWeightModel_thrPivEne", m_maxEdgeWeightModel_thrPivEne= 1000.);
  declareProperty ("m_maxEdgeWeightModel_thrHE"    , m_maxEdgeWeightModel_thrHE    = 200. );

}
    
StatusCode CalMSTClusteringTool::initialize()
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
    
    // check if jobOptions were passed
    log<<MSG::DEBUG<<"m_maxNumXtals                  set to "<<m_maxEdgeWeight  	     <<endreq;
    log<<MSG::DEBUG<<"m_maxEdgeWeightModel_thrLE     set to "<<m_maxEdgeWeightModel_thrLE    <<endreq;
    log<<MSG::DEBUG<<"m_maxEdgeWeightModel_thrPivEne set to "<<m_maxEdgeWeightModel_thrPivEne<<endreq;
    log<<MSG::DEBUG<<"m_maxEdgeWeightModel_thrHE     set to "<<m_maxEdgeWeightModel_thrHE    <<endreq;

    // Cluster filling utility
    // -- To be replaced with a generic version soon
    m_clusterInfo = new MomentsClusterInfo(m_calReconSvc);

    return StatusCode::SUCCESS ;
}

StatusCode CalMSTClusteringTool::findClusters(Event::CalClusterCol* calClusterCol)
{
    //Purpose and method:
    //
    //   This function performs the calorimeter clustering using MST's
    //   The main actions are:
    //      - calculate Mimimum Spanning Tree with all the xtals
    //                  (http://en.wikipedia.org/wiki/Minimum_spanning_tree)
    //                  using Prim's algorithm
    //      - fill UbreTree and UberEdged
    //      - remove all the edges whose weight is above threshold
    //                  and fill the list of Trees
    //      - store all calculated quantities in CalCluster objects
    // 
    // TDS input: CalXtalRecCol
    // TDS output: CalClustersCol

    // prepare the initital set of xtals
    //XtalDataList* xTalClus = new XtalDataList();

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

    // Make sure we have clean vectors to begin
    m_xTals.clear();
    m_xTals_setA.clear();
    m_uberTree.clear();
    m_clusterTree.clear();
    m_uberEdges.clear();

    // get list of xtals 
    for (Event::CalXtalRecCol::const_iterator it = m_calReconSvc->getXtalRecs()->begin() ; 
                it != m_calReconSvc->getXtalRecs()->end(); ++it )
    {
        // get pointer to the reconstructed data for given crystal
        Event::CalXtalRecData * recData = *it ;
       
	XtalDataListIterator xTalIter = m_xTals_setA.insert(m_xTals_setA.end(), recData);
    }

    // DEBUG Print
    // take a look at list of clusters.
    //std::cout << "RRRRRRRRRRRRRRRRRRR Total number of xtals: " <<m_xTals_setA.size() << std::endl;
    // End of DEBUG Print

    // case of 0 and 1 xtal will be handled separately...TBD
    if (m_xTals_setA.size()>1)
      {

	// Loop to fill the Uber MST
	
	// select a starting xtal
	m_xTals.push_back(m_xTals_setA.front());
	m_xTals_setA.pop_front();

	// DEBUG Print
	//std::cout << "RRRRRRRRRRRRRRRRRRR We start the loop with setA size = " 
	//	  << m_xTals_setA.size() << " and setB size = " << m_xTals.size() << std::endl;
	//std::cout << "WWWWWWWWWWW Filling the Uber Tree:" << std::endl;
	//std::cout << "Map1\tE1\tX1\tY1\tZ1\tMap2\tE2\tX2\tY2\tZ2\tW" << std::endl;
	// End of DEBUG Print

	// loop until there are unassociated xtals
	while(! m_xTals_setA.empty())
	  {
	    double minWeight = 100000000.; // init to a very long number
	    Event::CalXtalRecData * bestXtal1 = 0;
	    Event::CalXtalRecData * bestXtal2 = 0;
	    // loop over m_xTals
	    for (XtalDataList::iterator it=m_xTals.begin() ; it != m_xTals.end(); it++ )
	      {
		// then loop over m_xtals_setA
		for (XtalDataList::iterator itA=m_xTals_setA.begin() ; itA != m_xTals_setA.end(); itA++ )
		  {
		    // find the closest link between the two sets and the corresponding xtals
		    double currWeight = xtalsWeight( *it, *itA);
		    if (currWeight < minWeight)
		      {
			minWeight = currWeight;
			bestXtal1 = *it;
			bestXtal2 = *itA;
		      }
		  }
	      }
	      
	    // Fill the uber tree
	    m_uberEdges.push_back( MSTEdge(*bestXtal1,*bestXtal2, sqrt(minWeight)) ); 
	    m_uberTree.addEdge( m_uberEdges.back());

	    // DEBUG Print
	    //idents::CalXtalId xTalId1 = bestXtal1->getPackedId();
	    //idents::CalXtalId xTalId2 = bestXtal2->getPackedId();	    
	    //Point xTalPoint1 = bestXtal1->getPosition();
	    //Point xTalPoint2 = bestXtal2->getPosition();
	    ///* std::cout << "Map1\tE1\tX1\tY1\tZ1\tMap2\tE2\tX2\tY2\tZ2\tW" << std::endl; */
	    //std::cout <<  xTalId1.getPackedId() <<"\t" << bestXtal1->getEnergy()  <<"\t"<< xTalPoint1.x()  <<"\t"<< xTalPoint1.y()  <<"\t"<< xTalPoint1.z()  <<"\t";
	    //std::cout <<  xTalId2.getPackedId() <<"\t" << bestXtal2->getEnergy()  <<"\t"<< xTalPoint2.x()  <<"\t"<< xTalPoint2.y()  <<"\t"<< xTalPoint2.z()  <<"\t";
	    //std::cout << sqrt(minWeight) << std::endl;
	    // end of DEBUG Print

	    
	    // add best xtal2 to m_xTals and remove it from m_xTal_setA
	    m_xTals.push_back(bestXtal2);
    
	    // Now remove best xtal2 the crystal from m_xTals_setA 
	    // from CalSimpleClusteringTool::removeXTal ??? Why Tracy made this thing so complicated ?
	    XtalDataList::iterator xTalDataIter = std::find(m_xTals_setA.begin(), m_xTals_setA.end(), bestXtal2);
	    if (xTalDataIter != m_xTals_setA.end())
	      {
		m_xTals_setA.erase(xTalDataIter);
	      }
	    else
	      {
		int j = 0; 
	      }
	  }


	// calculate stats for uber tree
	m_uberTree.evalStats();
	m_maxEdgeWeight = getWeightThreshold(m_uberTree.getTotalEnergy());

	// DEBUG Print
	//std::cout << "WWWWWWWWWWWSSSSSSSSSSS Final Uber Tree size: " << m_uberTree.size() << std::endl;
	//std::cout << "WWWWWWWWWWWSSSSSSSSSSS Max weight for this event is: " << m_maxEdgeWeight << " for E= " << m_uberTree.getTotalEnergy() << std::endl;
	// end of DEBUG Print

	// Now we have the uber tree, need to loop over its edges 
	// and remove those above threshold;
	
	// create a first tree
	m_clusterTree.push_back(MSTTree());
	
	std::list<MSTEdge*> uberEdges = m_uberTree.getEdges(); // 
	int myEdgeCounter = 0;
	for (std::list<MSTEdge*>::iterator it=uberEdges.begin();  it != uberEdges.end(); it++ )
	{
	  MSTEdge* thisEdge = *it;
	  if (thisEdge->getWeight() > m_maxEdgeWeight) // time to split the tree.
	    {
	      // DEBUG Print
	      //idents::CalXtalId xTalId1 = thisEdge->getNode1()->getPackedId();
	      //idents::CalXtalId xTalId2 = thisEdge->getNode2()->getPackedId();
	      //std::cout << "WWWWWWWWWWWCCCCCCCC Found a large weight: " << thisEdge->getWeight() << std::endl;
	      //std::cout << "WWWWWWWWWWWCCCCCCCC Node1: Map " << xTalId1.getPackedId() << " E="<<  thisEdge->getNode1()->getEnergy() << std::endl;
	      //std::cout << "WWWWWWWWWWWCCCCCCCC Node2: Map " << xTalId2.getPackedId() << " E="<<  thisEdge->getNode2()->getEnergy() << std::endl;
	      // end of DEBUG Print

	      // Add node ONLY if splitting the very first edge.
	      if (myEdgeCounter == 0){m_clusterTree.back().addNode(thisEdge->getNode1());}

	      // a new tree is created
	      m_clusterTree.push_back(MSTTree());
	      m_clusterTree.back().addNode(thisEdge->getNode2()); 
	    }
	  else
	    {
	      m_clusterTree.back().addEdge(*thisEdge);
	    }
	  myEdgeCounter+=1;
	}


	// Sorting Trees - required since there is no intrinsic ordering in the MST algorithm
	// First eval some stats - energy and edges properties - fill members
	for (std::list<MSTTree>::iterator it= m_clusterTree.begin();  it !=  m_clusterTree.end(); it++ )
	  {
	    it->evalStats();	    
	  }
	// Second, real sorting.
	m_clusterTree.sort(compare_total_energy);

	// DEBUG Print
	// print all the nodes in the tree to check that is works.
	//int count = 0;
	//std::cout << "WWWWWWWWWWWCCCCCCCC Number of trees " << m_clusterTree.size() << std::endl;
	//for (std::list<MSTTree>::iterator it= m_clusterTree.begin();  it !=  m_clusterTree.end(); it++ )
	//  {
	//    std::cout << "WWWWWWWWWWWCCCCCCCC Tree number " << count << std::endl;
	//    count++;
	//    MSTTree thisTree = *it;
	//    std::cout << "WWWWWWWWWWWCCCCCCCC TotalEnergy " << thisTree.getTotalEnergy() << " MaxEnergy " << thisTree.getMaxEnergy() << std::endl;
	//    thisTree.printNodes();
	//  }
	// end of DEBUG Print
	
      

	// Now add the uber tree to the end of our list
	// But only do so if more than one cluster found
	if (m_clusterTree.size() > 1) m_clusterTree.push_back(m_uberTree);
    
	//-------------------------------------
	// Convert the results into CalClusters
	calClusterCol->clear() ;

	for (std::list<MSTTree>::iterator treeIter = m_clusterTree.begin(); treeIter != m_clusterTree.end(); treeIter++)
	  {
	    // create a temporary list of xtals
	    XtalDataList *xTalClus = new XtalDataList();
	    // fill the temporary list of xtals with the xtals in tree map.
	    std::map<int, Event::CalXtalRecData*> xtalMap = treeIter->getNodeMap();
	    std::map<int, Event::CalXtalRecData*>::const_iterator mapIter;
	    for ( mapIter=xtalMap.begin() ; mapIter != xtalMap.end(); mapIter++ )
	      {
		// get pointer to the reconstructed data for given crystal
		Event::CalXtalRecData * recData = (*mapIter).second ;
		xTalClus->push_back(recData);

		// end of DEBUG Print	
		//std::cout << "--------> FILLING CLUSTERS MAP: "<< (*mapIter).first << " => " << recData->getEnergy() << " MeV"<< std::endl;

	      }

	    // create and fill the cluster - from Tracy
	    Event::CalCluster* cluster = m_clusterInfo->fillClusterInfo(xTalClus);
	    std::string producerName("CalMSTClusteringTool/") ;
	    producerName += cluster->getProducerName() ;
	    cluster->setProducerName(producerName) ;
	    cluster->clearStatusBit(Event::CalCluster::ALLXTALS);
	    
	    // Set the CalMSTreeParams for the cluster and associated status bit    
	    // make sure first that statistics are available
	    if(treeIter->getStatsBit()==false) treeIter->evalStats();
	    Event::CalMSTreeParams mstreeparams(treeIter->getTotalEnergy(),
 		  treeIter->getMaxEnergy(),	treeIter->getNumEdges(),		     
 		  treeIter->getMinEdgeLength(), treeIter->getMaxEdgeLength(),		     
 		  treeIter->getMeanEdgeLength(),treeIter->getMeanEdgeLengthTrunc(),	     
 		  treeIter->getRmsEdgeLength(), treeIter->getRmsEdgeLengthTrunc());
		  
	    cluster->setMSTreeParams(mstreeparams);
	    cluster->setStatusBit(Event::CalCluster::MSTTREE);
	    // Raw debugging -- works
	    //std::cout<<"CalMSTreeParams "<<cluster->getMSTreeParams()<<std::endl;
	    	    
	    // Add cluster into the collection
	    calClusterCol->push_back(cluster);
	    // Do I need to delete mstreeparams ? -- Johan

	    // Loop through the xtals to make the relational table (hmmm... this could be done better...)  - from Tracy
	    for(XtalDataListIterator xTalIter = xTalClus->begin(); xTalIter != xTalClus->end(); xTalIter++)
	      {
		Event::CalXtalRecData*   xTal         = *xTalIter;
		Event::CalClusterHitRel* xTal2ClusRel = new Event::CalClusterHitRel(xTal,cluster);
		
		m_xTal2ClusTab->addRelation(xTal2ClusRel);
	      }
	    
	    delete xTalClus;
	  }
	
      } // end of case (m_xTals_setA.size()>1)
    else // special handling for m_xTals_setA.size() = 1 or 0
      {
	// what follows comes form calSingleClusteringTool
	calClusterCol->clear() ;
	
	// Get the cluster instance - only one obviously
	Event::CalCluster* cluster = m_clusterInfo->fillClusterInfo(&m_xTals_setA);

	std::string producerName("CalMSTClusteringTool/") ;
	producerName += cluster->getProducerName() ;
	cluster->setProducerName(producerName) ;
	cluster->setStatusBit(Event::CalCluster::ALLXTALS);
	
	// Attach a dummy mstree params container - but do not set the status bit
	Event::CalMSTreeParams mstreeparams(-1,-1,-1,-1,-1,-1,-1,-1,-1);
	cluster->setMSTreeParams(mstreeparams);
	 
	calClusterCol->push_back(cluster);

	// Loop through again to make the relations
	for (Event::CalXtalRecCol::const_iterator it = m_calReconSvc->getXtalRecs()->begin() ; 
	     it != m_calReconSvc->getXtalRecs()->end(); ++it )
	  {
	    // get pointer to the reconstructed data for given crystal
	    Event::CalXtalRecData * recData = *it ;

	    // Even though only one cluster, we still need to make the relations!
	    Event::CalClusterHitRel* xTal2ClusRel = new Event::CalClusterHitRel(recData,cluster);

	    xTal2ClusTabList->push_back(xTal2ClusRel);
	  }
	
      } // end of clustering
 

    // And delete the table (remembering that the list is in the TDS)
    delete m_xTal2ClusTab;
 

    return StatusCode::SUCCESS ;
}


// weight calculation between 2 xtals
double CalMSTClusteringTool::xtalsWeight(Event::CalXtalRecData* xTal1, Event::CalXtalRecData* xTal2 )
{
  // calculate the weights of the line that connects two xtals
  // now is the (square of the) euclidean distance, but is can be a more complex funcion
  // keep in mind that this function is called O(N^2) times per events (N= number od xtals)
  // so it must be as fast as possible.

  //Compute distance to this xTal

  // is this faster? 
  //Vector distVec = xTal1->getPosition() - xTal2->getPosition();
  //double dist2 = distVec.square();

  Point xTalPoint1 = xTal1->getPosition();
  Point xTalPoint2 = xTal2->getPosition();

  double dist2 = (xTalPoint1.x() - xTalPoint2.x())*(xTalPoint1.x() - xTalPoint2.x()) +
    (xTalPoint1.y() - xTalPoint2.y())*(xTalPoint1.y() - xTalPoint2.y()) +
    (xTalPoint1.z() - xTalPoint2.z())*(xTalPoint1.z() - xTalPoint2.z());

  return dist2 ;
}


double CalMSTClusteringTool::getWeightThreshold(double energy)
{
  // Implement a max weight per events that is energy dependent, defaults:
  //  m_maxEdgeWeightModel_thrLE    =500. 
  //  m_maxEdgeWeightModel_thrPivEne=1000.
  //  m_maxEdgeWeightModel_thrHE    =200. 

  if (energy>m_maxEdgeWeightModel_thrPivEne)
      return m_maxEdgeWeightModel_thrHE;
  else
      return ( m_maxEdgeWeightModel_thrLE - 150.0*(log10(energy) - 1));
  
}

