
#include <reconRootData/CalCluster.h>
#include <commonRootData/RootDataUtil.h>
#include "Riostream.h"

ClassImp(CalCluster)

void CalCluster::init
 ( const std::vector<CalClusterLayerData> & layers,
   const CalMSTreeParams & treeParams,
   const CalFitParams & fitParams, 
   const CalMomParams & momParams,
   const CalClassParams & classParams,
   Int_t numSaturatedXtals, Int_t numTruncXtals, UInt_t statusBits )
 {
  std::vector<CalClusterLayerData>::const_iterator layer ;
  int i ;
  for ( i=0, layer=layers.begin() ; layer != layers.end() ; ++i, ++layer )
    m_layers[i] = (*layer) ;
  m_mstParams    = treeParams;
  m_fitParams    = fitParams;
  m_momParams    = momParams ;
  m_classParams  = classParams;
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
  m_classParams.Clear();
  m_numTruncXtals = 0 ;
  m_statusBits = 0 ;
 }

CalCluster::CalCluster
 ( const std::vector<CalClusterLayerData> & layers,
   const CalMSTreeParams & treeParams,
   const CalFitParams& fitParams,
   const CalMomParams & momParams,
   const CalClassParams & classParams,
   Int_t numSaturatedXtals, Int_t numTruncXtals, UInt_t statusBits )
 {
   init(layers,treeParams,fitParams,momParams,classParams,
	numSaturatedXtals,numTruncXtals,statusBits) ;
 }

Double_t CalCluster::getClassProb(const std::string& className) const
{
  return m_classParams.getProb(className);
}

Double_t CalCluster::getGamProb() const
{
  return getClassProb("gam");
}

void CalCluster::Print( Option_t * ) const
 {
   std::cout << 
     "Cal cluster status bits: " << m_statusBits << "\n" <<
     "Number of saturated xtals: " << m_numSaturatedXtals << "\n" <<
     "Truncated number of xtals: " << m_numTruncXtals << "\n" <<
     "----------------------------------------------------\n" <<
     "---------- Output from the MST clustering ----------\n" << 
     "----------------------------------------------------\n";
   m_mstParams.Print() ;
   std::cout << 
     "----------------------------------------------------\n" <<
     "----------- Output from the fitting tool -----------\n" << 
     "----------------------------------------------------\n";
   m_fitParams.Print() ;
   std::cout << 
     "----------------------------------------------------\n" <<
     "--------- Output from the moments analysis ---------\n" << 
     "----------------------------------------------------\n";
   m_momParams.Print() ;
   std::cout << 
     "------------------------------------------------------" << std::endl;
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
    CalClassParams c;
    c.Fake(ievent,rank,randNum) ;
    CalClusterLayerData layer ;
    std::vector<CalClusterLayerData> clusLayerData ;
    UInt_t iclusLayer;
    for ( iclusLayer = 0; iclusLayer<ROOT_NUMCALLAYERS ; ++iclusLayer ) {
        layer.Fake(ievent,rank*ROOT_NUMCALLAYERS+iclusLayer,randNum) ;
        clusLayerData.push_back(layer) ;
    }

    init(clusLayerData,t,f,p,c,4,4,5) ;
 }

Bool_t CalCluster::CompareInRange( const CalCluster & c, const std::string & name ) const {

    bool result = true ;

    result = getMSTreeParams().CompareInRange(c.getMSTreeParams()) && result ;
    result = getFitParams().CompareInRange(c.getFitParams()) && result ;
    result = getMomParams().CompareInRange(c.getMomParams()) && result ;
    result = getClassParams().CompareInRange(c.getClassParams()) && result ;

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

