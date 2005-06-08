#include "reconRootData/CalCluster.h"
#include <iostream>

ClassImp(CalCluster)

void CalCluster::init
 ( const std::vector<CalClusterLayerData> & layers,
   const CalParams & params,
   Double_t rmsLong, Double_t rmsTrans )
 {
  std::vector<CalClusterLayerData>::const_iterator layer ;
  int i ;
  for ( i=0, layer=layers.begin() ; layer != layers.end() ; ++i, ++layer )
    m_layers[i] = (*layer) ;
  m_params = params ;
  m_rmsLong = rmsLong ;
  m_rmsTrans = rmsTrans ;
 }

void CalCluster::Clear( Option_t * )
 {
  int i ;
  for ( i=0 ; i<ROOT_NUMCALLAYERS ; ++i )
    m_layers[i].Clear() ;
  m_params.Clear() ;
  m_rmsLong = 0.0 ;
  m_rmsTrans = 0.0 ;
 }

CalCluster::CalCluster
 ( const std::vector<CalClusterLayerData> & layers,
   const CalParams & params,
   Double_t rmsLong, Double_t rmsTrans )
 { init(layers,params,rmsLong,rmsTrans) ; }

void CalCluster::Print( Option_t * ) const
 {
  m_params.Print() ;
  using namespace std ;
  cout<<"RMS Long: "<<m_rmsLong<<"  RMS Trans: "<<m_rmsTrans<<endl ;
  cout<<"Number of layers: "<<ROOT_NUMCALLAYERS<<endl ;
 }

