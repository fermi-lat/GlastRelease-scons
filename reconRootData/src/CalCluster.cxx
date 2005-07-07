
#include <reconRootData/CalCluster.h>
#include <commonRootData/RootDataUtil.h>
#include <iostream>

ClassImp(CalCluster)

void CalCluster::init
 ( const std::vector<CalClusterLayerData> & layers, const CalParams & params,
   Double_t rmsLong, Double_t rmsLongAsym, Double_t rmsTrans,
   Int_t numTruncXtals, UInt_t statusBits )
 {
  std::vector<CalClusterLayerData>::const_iterator layer ;
  int i ;
  for ( i=0, layer=layers.begin() ; layer != layers.end() ; ++i, ++layer )
    m_layers[i] = (*layer) ;
  m_params = params ;
  m_rmsLong = rmsLong ;
  m_rmsLongAsym = rmsLongAsym ;
  m_rmsTrans = rmsTrans ;
  m_numTruncXtals = numTruncXtals ;
  m_statusBits = statusBits ;
 }

void CalCluster::Clear( Option_t * )
 {
  int i ;
  for ( i=0 ; i<ROOT_NUMCALLAYERS ; ++i )
    m_layers[i].Clear() ;
  m_params.Clear() ;
  m_rmsLong = 0.0 ;
  m_rmsLongAsym = 0.0 ;
  m_rmsTrans = 0.0 ;
  m_numTruncXtals = 0 ;
  m_statusBits = 0 ;
 }

CalCluster::CalCluster
 ( const std::vector<CalClusterLayerData> & layers, const CalParams & params,
   Double_t rmsLong, Double_t rmsLongAsym, Double_t rmsTrans,
   Int_t numTruncXtals, UInt_t statusBits )
 { init(layers,params,rmsLong,rmsLongAsym,rmsTrans,numTruncXtals,statusBits) ; }

// dummy data, just for tests
void CalCluster::Fake( UInt_t rank, Float_t randNum )
 {
    CalParams p ;
    p.Fake(rank,randNum) ;
    CalClusterLayerData layer ;
    std::vector<CalClusterLayerData> clusLayerData ;
    UInt_t iclusLayer;
    for ( iclusLayer = 0; iclusLayer<ROOT_NUMCALLAYERS ; ++iclusLayer ) {
        layer.Fake(rank*ROOT_NUMCALLAYERS+iclusLayer,randNum) ;
        clusLayerData.push_back(layer) ;
    }
    init(clusLayerData,p,1.0,2.0,3.0,4,5) ;
 }

void CalCluster::Print( Option_t * ) const
 {
  m_params.Print() ;
  std::cout
    <<"No.Trunc Xtals "<<m_numTruncXtals<<"\n"
    <<"RMS Long: "<<m_rmsLong<<"RMS Long Asym: "<<m_rmsLongAsym<<"  RMS Trans: "<<m_rmsTrans<<"\n"
    <<"Number of layers: "<<ROOT_NUMCALLAYERS<<std::endl ;
 }

Bool_t CalCluster::CompareInRange( const CalCluster & c ) const {

    bool result = true ;

    result = result && getParams().CompareInRange(c.getParams()) ;
    int i ;
    for ( i=0 ; i<ROOT_NUMCALLAYERS ; ++i ) {
        result = result && getLayer(i).CompareInRange(c.getLayer(i)) ;
    }

    result = result && rootdatautil::CompareInRange(getRmsLong(),c.getRmsLong(),"RmsLong") ;
    result = result && rootdatautil::CompareInRange(getRmsLongAsym(),c.getRmsLongAsym(),"RmsLongAsym") ;
    result = result && rootdatautil::CompareInRange(getRmsTrans(),c.getRmsTrans(),"RmsTrans") ;
    result = result && rootdatautil::CompareInRange(getNumTruncXtals(),c.getNumTruncXtals(),"NumTruncXtals") ;
    result = result && rootdatautil::CompareInRange(getStatusBits(),c.getStatusBits(),"StatusBits") ;

    if (!result) {
        std::cout<<"Comparison ERROR for "<<ClassName()<<std::endl ;
    }
    return result ;

}

