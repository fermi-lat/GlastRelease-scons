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

#include <CalRecon/ICalReconSvc.h>
#include <CalRecon/ICalClusteringTool.h>
#include "StdClusterInfo.h"
#include "MomentsClusterInfo.h"

#include "TMath.h"

// A useful typedef 
typedef  XtalDataList::iterator XtalDataListIterator;

// Define the map between the layer/row index and the crystals
// typedef std::map <int, XtalDataList> LyrRow2XtalDataMap; // Tracy stuff

//
// Define a class for MSTEdge
// An Edge is a weight and two xtals (connected by the line with a weight)
//

class MSTEdge
{

public:
  MSTEdge( )
    :node1(0),node2(0),weight(0.){;};
  MSTEdge(const Event::CalXtalRecData& cry1, const Event::CalXtalRecData& cry2, double w)
    :node1(&cry1),node2(&cry2),weight(w){;}
  MSTEdge(const MSTEdge& other)
    :node1(other.node1),node2(other.node2),weight(other.weight){;}
  ~MSTEdge( )  { };
  
  void setEdge(const Event::CalXtalRecData*, const Event::CalXtalRecData*, const double);
  const Event::CalXtalRecData* getNode1(){return node1;};
  const Event::CalXtalRecData* getNode2(){return node2;};
  double getWeight(){return weight;};
  // stupid-useless function that I use for debug
  double getSumEnergy() { return (node1->getEnergy() + node2->getEnergy()) ;} 
  void printSumEnergy() {std::cout << " Ene1 = " << node1->getEnergy() << " Ene2 = " << node2->getEnergy() << std::endl;};

private:
  const Event::CalXtalRecData* node1;
  const Event::CalXtalRecData* node2;
  double weight;

};

void MSTEdge::setEdge(const Event::CalXtalRecData* _node1, const  Event::CalXtalRecData* _node2, const double _weight)
{
  node1 = _node1;
  node2 = _node2;
  weight = _weight;
}

// MST Tree: a list of edges + something that may be useful later
// should I use pointers?
// should I define a destructor and explicitely kill all the objects?
class MSTTree
{

public:
  MSTTree( ) ;
  ~MSTTree( )  {};
  
  //void addEdge(MSTEdge &_edg) {edges.push_back(&_edg);}; // this works
  void addEdge(MSTEdge &);
  void addNode(const Event::CalXtalRecData* _node) ; // implemented using a map
  int size () {return (edges.size());}
  void clear() {edges.clear(); nodemap.clear();};
  std::list<MSTEdge*> getEdges() {return edges;}; 
  std::map<int, const  Event::CalXtalRecData*> getNodeMap() {return nodemap;};
  void printNodes();
  
private:
  int  getXtalUniqueId(const Event::CalXtalRecData *);
  std::list<MSTEdge*> edges;
  double totalEnergy;
  //const Event::CalXtalRecData* node;
  std::map<int, const  Event::CalXtalRecData*> nodemap;
  
};

MSTTree::MSTTree () 
{
  
  edges.clear(); // make sure with know the initial state.
  totalEnergy = 0.;
  //node = NULL;
  nodemap.clear();
}

void MSTTree::addEdge(MSTEdge &_edg) 
{
  edges.push_back(&_edg);
  // add nodes
  //addNode(_edg.getNode1()); // Segmentation fault?
  //addNode(_edg.getNode2()); //
  
  addNode(edges.back()->getNode1()); 
  addNode(edges.back()->getNode1()); 
}

void MSTTree::addNode(const Event::CalXtalRecData *_node) 
{
  //std::cout << "xtal LyrRow " << getXtalUniqueId(_node) << std::endl;
  nodemap[getXtalUniqueId(_node)] = _node;
}

void MSTTree::printNodes() 
{
  // example of show content:
  std::map<int, const  Event::CalXtalRecData*>::iterator it;
  for ( it=nodemap.begin() ; it != nodemap.end(); it++ )
    std::cout << "--------> MAP: "<< (*it).first << " => " << (*it).second->getEnergy() << " MeV"<< std::endl;

}


// a unique Id to avodi duplication in the nodemap
int  MSTTree::getXtalUniqueId(const Event::CalXtalRecData * xTal)
{
  idents::CalXtalId xTalId = xTal->getPackedId();
        
    //    int curLayer  = xTalId.getLayer();
    //int curTower  = xTalId.getTower();
    //int row       = xTalId.isX() ? curTower % 4 : curTower / 4;
    // int layerRow  = 4 * curLayer + row;
   
  //int layerRow  = *xTal+1000;
  //return layerRow;
  return xTalId.getPackedId();
}


// end of my test classes

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
  
  // calculate thw weight between two xtals
  double xtalsWeight(Event::CalXtalRecData* xTal1, Event::CalXtalRecData* xTal2 );
  // uses the Xtal identifier to build a layer/row index
  int  getXtalLayerRow(Event::CalXtalRecData* xTal);
  
  //! Service for basic Cal info
  ICalReconSvc*      m_calReconSvc;
  
  //! Event Service member directly useable by concrete classes.
  IDataProviderSvc*  m_dataSvc;
  
  //! Utility for filling clusters
  ICalClusterFiller* m_clusterInfo;

  // Keep trackk of crystals locally
  XtalDataList             m_xTals;
  XtalDataList             m_xTals_setA;

  // MSTtree as list of edges
  //std::list<MSTEdge> m_uberTree;
  MSTTree m_uberTree;
  std::list<MSTEdge>  m_uberEdges;

  // Clusters are based on a list of trees
  std::list<MSTTree> m_clusterTree;

  // Parameter to separate trees
  double       m_maxEdgeWeight;
} ;

DECLARE_TOOL_FACTORY(CalMSTClusteringTool) ;

CalMSTClusteringTool::CalMSTClusteringTool(const std::string & type, 
                                                 const std::string & name, 
                                                 const IInterface* parent)
  : AlgTool(type,name,parent),
    m_maxEdgeWeight(300.)
{ 

  m_xTals.clear();
  m_xTals_setA.clear();
  m_clusterTree.clear();
  m_uberEdges.clear();
  declareInterface<ICalClusteringTool>(this) ; 
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

    // Cluster filling utility
    // -- To be replaced with a generic version soon
    m_clusterInfo = new MomentsClusterInfo(m_calReconSvc);

    return StatusCode::SUCCESS ;
}

StatusCode CalMSTClusteringTool::findClusters(Event::CalClusterCol* calClusterCol)
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
    XtalDataList* xTalClus = new XtalDataList();

    // Create a Xtal to Cluster relations list
    Event::CalClusterHitTabList* xTal2ClusTabList = new Event::CalClusterHitTabList();
    xTal2ClusTabList->clear();

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
    xTalClus->clear();
    for (Event::CalXtalRecCol::const_iterator it = m_calReconSvc->getXtalRecs()->begin() ; 
                it != m_calReconSvc->getXtalRecs()->end(); ++it )
    {
        // get pointer to the reconstructed data for given crystal
        Event::CalXtalRecData * recData = *it ;
        xTalClus->push_back(recData) ;

	XtalDataListIterator xTalIter = m_xTals_setA.insert(m_xTals_setA.end(), recData);
    }

    // take a look at list of clusters.
    std::cout << "RRRRRRRRRRRRRRRRRRR Total number of xtals: " <<m_xTals_setA.size() << std::endl;

    // case of 1 xtal will be handled separately...TBD
    if (m_xTals_setA.size()>1 && m_xTals_setA.size()<100)
      {

	// Loop to fill the Uber MST
	
	// select a starting xtal
	m_xTals.push_back(m_xTals_setA.front());
	m_xTals_setA.pop_front();
	std::cout << "RRRRRRRRRRRRRRRRRRR We start the loop with setA size = " 
		  << m_xTals_setA.size() << " and setB size = " << m_xTals.size() << std::endl;
	// loop until there are unassociated xtals
	while(! m_xTals_setA.empty())
	  {
	    double minWeight = 100000000000; // init to a very long number
	    Event::CalXtalRecData * bestXtal1;
	    Event::CalXtalRecData * bestXtal2;
	    // loop over m_xTals
	    for (XtalDataList::iterator it=m_xTals.begin() ; it != m_xTals.end(); it++ )
	      {
		// then loop over m_xtals_setA
		for (XtalDataList::iterator itA=m_xTals_setA.begin() ; itA != m_xTals_setA.end(); itA++ )
		  {
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
	    m_uberEdges.back().printSumEnergy();
	    m_uberTree.addEdge( m_uberEdges.back());

	    m_uberTree.printNodes();
	    std::cout << "WWWWWWWWWWW Filling the Uber Tree with weight: " << sqrt(minWeight) << std::endl;

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
	    // just a test output to check we do not have an infinite loop. 
	    //std::cout << "WWWWWWWWWWW Check that setA size is decreasing: " <<  m_xTals_setA.size() << " and setB size = " << m_xTals.size() << std::endl;
	  }
	std::cout << "WWWWWWWWWWWSSSSSSSSSSS Final Uber Tree size: " << m_uberTree.size() << std::endl;

	// Now we have the uber tree, need to loop over its edges 
	// and remove those above threshold;
	
	m_clusterTree.push_back(MSTTree());
	

	std::list<MSTEdge*> uberEdges = m_uberTree.getEdges(); // 
	for (std::list<MSTEdge*>::iterator it=uberEdges.begin();  it != uberEdges.end(); it++ )
	{
	  MSTEdge* thisEdge = *it;
	  if (thisEdge->getWeight() > m_maxEdgeWeight)
	    {
	      std::cout << "WWWWWWWWWWWCCCCCCCC Found a large weight: " << thisEdge->getWeight() << std::endl;
	      m_clusterTree.back().addNode(thisEdge->getNode1());
	      // a new tree is created
	      m_clusterTree.push_back(MSTTree());
	      m_clusterTree.back().addNode(thisEdge->getNode2()); 
	    }
	  else
	    {
	      m_clusterTree.back().addEdge(*thisEdge);
	    }
	  
	}

	// print all the nodes in the tree to check that is works.
	int count = 0;
	std::cout << "WWWWWWWWWWWCCCCCCCC Number of trees " << m_clusterTree.size() << std::endl;
	for (std::list<MSTTree>::iterator it= m_clusterTree.begin();  it !=  m_clusterTree.end(); it++ )
	  {
	    std::cout << "WWWWWWWWWWWCCCCCCCC Tree number " << count << std::endl;
	    count++;
	    MSTTree thisTree = *it;
	    thisTree.printNodes();
	  }
	
      };



    calClusterCol->clear() ;

    // Get the cluster instance
    Event::CalCluster* cluster = m_clusterInfo->fillClusterInfo(xTalClus);

    std::string producerName("CalMSTClusteringTool/") ;
    producerName += cluster->getProducerName() ;
    cluster->setProducerName(producerName) ;
    cluster->setStatusBit(Event::CalCluster::ALLXTALS); 
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

    delete xTalClus;
  
    return StatusCode::SUCCESS ;
}

// from Tracy simple clustering tool
int  CalMSTClusteringTool::getXtalLayerRow(Event::CalXtalRecData* xTal)
{
    idents::CalXtalId xTalId = xTal->getPackedId();
        
    int curLayer  = xTalId.getLayer();
    int curTower  = xTalId.getTower();
    int row       = xTalId.isX() ? curTower % 4 : curTower / 4;
    int layerRow  = 4 * curLayer + row;

    return layerRow;
}

double CalMSTClusteringTool::xtalsWeight(Event::CalXtalRecData* xTal1, Event::CalXtalRecData* xTal2 )
{
  // calculate the weights of the line that connects two xtals
  // now is the (square of the) euclidean distance, but is can be a more complex funcion
  // keep in mind that this function is called O(N^2) times per events (N= number od xtals)
  // so it must be as fast as possible.

  //Compute distance to this xTal
  Vector distVec = xTal1->getPosition() - xTal2->getPosition();
  double dist2 = distVec.square();
  
  return dist2 ;
}
