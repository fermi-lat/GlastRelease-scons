
#include <reconRootData/CalCluster.h>
#include <commonRootData/RootDataUtil.h>
#include "Riostream.h"

ClassImp(CalCluster)

void CalCluster::init
 ( const std::vector<CalClusterLayerData> & layers,
   const CalMSTreeParams & treeParams,
   const CalFitParams & fitParams, 
   const CalMomParams & momParams,
   const std::map <std::string, double>& classprob,
   Int_t numSaturatedXtals, Int_t numTruncXtals, UInt_t statusBits )
 {
  std::vector<CalClusterLayerData>::const_iterator layer ;
  int i ;
  for ( i=0, layer=layers.begin() ; layer != layers.end() ; ++i, ++layer )
    m_layers[i] = (*layer) ;
  m_mstParams    = treeParams;
  m_fitParams    = fitParams;
  m_momParams    = momParams ;
  m_classesProb  = classprob;
  m_numSaturatedXtals = numSaturatedXtals;
  m_numTruncXtals = numTruncXtals ;
  m_statusBits = statusBits ;
 }

void CalCluster::Clear( Option_t * )
 {
  int i ;
  for ( i=0 ; i<ROOT_NUMCALLAYERS ; ++i )
    m_layers[i].Clear() ;
  m_mstParams.Clear();
  m_fitParams.Clear();
  m_momParams.Clear() ;
  m_classesProb.clear();
  m_numTruncXtals = 0 ;
  m_statusBits = 0 ;
 }

CalCluster::CalCluster
 ( const std::vector<CalClusterLayerData> & layers,
   const CalMSTreeParams & treeParams,
   const CalFitParams& fitParams, 
   const CalMomParams & momParams,
   const std::map <std::string, double>& classprob,
   Int_t numSaturatedXtals, Int_t numTruncXtals, UInt_t statusBits )
 {
   init(layers,treeParams, fitParams, momParams, classprob,
	numSaturatedXtals,numTruncXtals,statusBits) ;
 }

/// Access any classification probabilities
Double_t CalCluster::getTopologyProb(std::string top) const
{
    if(m_classesProb.count(top))
      return m_classesProb.find(top)->second;
    else
      return -1;
}   


void CalCluster::Print( Option_t * ) const
 {
  m_momParams.Print() ;
  m_fitParams.Print() ;
  m_mstParams.Print() ;
  std::cout
    <<"Gam prob " << getGamProb() << "\n"
    <<"No. Saturated Xtals " << m_numSaturatedXtals << "\n"
    <<"No.Trunc Xtals "<<m_numTruncXtals<<"\n"
    <<"Number of layers: "<<ROOT_NUMCALLAYERS<<std::endl ;
 }

// dummy data, just for tests
void CalCluster::Fake( Int_t ievent, UInt_t rank, Float_t randNum )
 {
    CalMSTreeParams t;
    t.Fake(ievent,rank,randNum);
    CalFitParams f;
    f.Fake(ievent,rank,randNum);
    CalMomParams p ;
    p.Fake(ievent,rank,randNum) ;
    CalClusterLayerData layer ;
    std::vector<CalClusterLayerData> clusLayerData ;
    UInt_t iclusLayer;
    for ( iclusLayer = 0; iclusLayer<ROOT_NUMCALLAYERS ; ++iclusLayer ) {
        layer.Fake(ievent,rank*ROOT_NUMCALLAYERS+iclusLayer,randNum) ;
        clusLayerData.push_back(layer) ;
    }
    std::map <std::string, double> m;
    m["gam"]=0.5;
    
    init(clusLayerData,t,f,p,m,4,4,5) ;
 }

Bool_t CalCluster::CompareInRange( const CalCluster & c, const std::string & name ) const {

    bool result = true ;

    result = getMomParams().CompareInRange(c.getMomParams()) && result ;
    result = getFitParams().CompareInRange(c.getFitParams()) && result ;
    result = getMSTreeParams().CompareInRange(c.getMSTreeParams()) && result ;

    int i ;
    for ( i=0 ; i<ROOT_NUMCALLAYERS ; ++i ) {
        result = getLayer(i).CompareInRange(c.getLayer(i)) && result ;
    }

    result = rootdatautil::CompareInRange(getGamProb(),c.getGamProb(),"GamProb") && result ;    
    result = rootdatautil::CompareInRange(getNumTruncXtals(),c.getNumTruncXtals(),"NumTruncXtals") && result ;
    result = rootdatautil::CompareInRange(getStatusBits(),c.getStatusBits(),"StatusBits") && result ;

    if (!result) {
        if ( name == "" ) {
            std::cout<<"Comparison ERROR for "<<ClassName()<<std::endl ;
        }
        else {
            std::cout<<"Comparison ERROR for "<<name<<std::endl ;
        }
    }
    return result ;

}

