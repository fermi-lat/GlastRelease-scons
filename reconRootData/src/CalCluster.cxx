
#include <reconRootData/CalCluster.h>
#include <commonRootData/RootDataUtil.h>
#include "Riostream.h"

ClassImp(CalCluster)

void CalCluster::init
 ( const std::vector<CalClusterLayerData> & layers, const CalFitParams & fitParams, 
   const CalParams & params,
   Double_t rmsLong, Double_t rmsLongAsym, Double_t rmsTrans, Double_t longSkew,
   Int_t numSaturatedXtals, Int_t numTruncXtals, UInt_t statusBits )
 {
  std::vector<CalClusterLayerData>::const_iterator layer ;
  int i ;
  for ( i=0, layer=layers.begin() ; layer != layers.end() ; ++i, ++layer )
    m_layers[i] = (*layer) ;
  m_fitParams = fitParams;
  m_params = params ;
  m_rmsLong = rmsLong ;
  m_rmsLongAsym = rmsLongAsym ;
  m_rmsTrans = rmsTrans ;
  m_skewnessLong = longSkew ;
  m_numSaturatedXtals = numSaturatedXtals;
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
  m_skewnessLong = -9999. ;
  m_numTruncXtals = 0 ;
  m_statusBits = 0 ;
 }

CalCluster::CalCluster
 ( const std::vector<CalClusterLayerData> & layers, const CalFitParams& fitParams, 
   const CalParams & params,
   Double_t rmsLong, Double_t rmsLongAsym, Double_t rmsTrans, Double_t longSkew,
   Int_t numSaturatedXtals, Int_t numTruncXtals, UInt_t statusBits )
{ init(layers,fitParams, params,rmsLong,rmsLongAsym,rmsTrans,longSkew,numSaturatedXtals,numTruncXtals,statusBits) ; }

void CalCluster::Print( Option_t * ) const
 {
  m_params.Print() ;
  std::cout
    <<"No. Saturated Xtals " << m_numSaturatedXtals << "\n"
    <<"No.Trunc Xtals "<<m_numTruncXtals<<"\n"
    <<"RMS Long: "<<m_rmsLong<<"RMS Long Asym: "<<m_rmsLongAsym<<"  RMS Trans: "<<m_rmsTrans<<"\n"
    <<"Number of layers: "<<ROOT_NUMCALLAYERS<<std::endl ;
 }

// dummy data, just for tests
void CalCluster::Fake( Int_t ievent, UInt_t rank, Float_t randNum )
 {
    CalFitParams f;
    f.Fake(ievent,rank,randNum);
    CalParams p ;
    p.Fake(ievent,rank,randNum) ;
    CalClusterLayerData layer ;
    std::vector<CalClusterLayerData> clusLayerData ;
    UInt_t iclusLayer;
    for ( iclusLayer = 0; iclusLayer<ROOT_NUMCALLAYERS ; ++iclusLayer ) {
        layer.Fake(ievent,rank*ROOT_NUMCALLAYERS+iclusLayer,randNum) ;
        clusLayerData.push_back(layer) ;
    }
    init(clusLayerData,f,p,1.0,2.0,3.0,3.5,4,4,5) ;
 }

Bool_t CalCluster::CompareInRange( const CalCluster & c, const std::string & name ) const {

    bool result = true ;

    result = getParams().CompareInRange(c.getParams()) && result ;
    int i ;
    for ( i=0 ; i<ROOT_NUMCALLAYERS ; ++i ) {
        result = getLayer(i).CompareInRange(c.getLayer(i)) && result ;
    }

    result = rootdatautil::CompareInRange(getRmsLong(),c.getRmsLong(),"RmsLong") && result ;
    result = rootdatautil::CompareInRange(getRmsLongAsym(),c.getRmsLongAsym(),"RmsLongAsym") && result ;
    result = rootdatautil::CompareInRange(getRmsTrans(),c.getRmsTrans(),"RmsTrans") && result ;
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

