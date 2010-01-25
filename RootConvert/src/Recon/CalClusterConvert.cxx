
#include <RootConvert/Recon/CalClusterConvert.h>
#include <RootConvert/Recon/CalFitParamsConvert.h>
#include <RootConvert/Recon/CalParamsConvert.h>
#include <RootConvert/Utilities/Toolkit.h>

namespace RootPersistence {

void convert( const Event::CalCluster & tdsCluster, CalCluster & rootCluster )
 {
    std::vector<CalClusterLayerData> rootLayers ;
    std::vector<Event::CalClusterLayerData>::const_iterator tdsLayer ;
    if (tdsCluster.size()!=NUMCALLAYERS)
     {
     }
    for ( tdsLayer=tdsCluster.begin() ;
          tdsLayer!=tdsCluster.end() ;
          ++tdsLayer )
    {
        CalClusterLayerData layerRoot
          ( tdsLayer->getEnergy(),
            convert(tdsLayer->getPosition()),
            convert(tdsLayer->getRmsSpread())) ;
        rootLayers.push_back(layerRoot) ;
    }
      
    CalFitParams rootFitParams ;
    convert(tdsCluster.getFitParams(),rootFitParams) ;
      
    CalParams rootParams ;
    convert(tdsCluster.getCalParams(),rootParams) ;

    rootCluster.init
      ( rootLayers, rootFitParams, rootParams,
        tdsCluster.getRmsLong(),
        tdsCluster.getRmsLongAsym(),
        tdsCluster.getRmsTrans(),
        (Int_t)tdsCluster.getNumSaturatedXtals(),
        (Int_t)tdsCluster.getNumTruncXtals(),
        tdsCluster.getStatusBits() ) ;
 }
 
void convert( const CalCluster & rootCluster, Event::CalCluster & tdsCluster )
{
    Event::CalFitParams tdsFitParams;
    convert(rootCluster.getFitParams(),tdsFitParams) ;    
    Event::CalParams tdsParams ;
    convert(rootCluster.getParams(),tdsParams) ;    

    tdsCluster.initialize
      ( tdsFitParams,
        tdsParams,
        rootCluster.getRmsLong(),
        rootCluster.getRmsTrans(),
        rootCluster.getRmsLongAsym(),
        rootCluster.getNumSaturatedXtals(),
        rootCluster.getNumTruncXtals() ) ;
    
    UInt_t rootStatusBits = rootCluster.getStatusBits() ;
    tdsCluster.setStatusBits(rootStatusBits) ;

    tdsCluster.clear() ;
    int i ;
    for ( i=0 ; i<NUMCALLAYERS ; ++i )
    {
        const CalClusterLayerData & rootLayer = rootCluster.getLayer(i) ;
        Event::CalClusterLayerData tdsLayer
          ( rootLayer.getEnergy(),
            convertPoint(rootLayer.getPosition()),
            convertVector(rootLayer.getRmsSpread()) ) ;
        tdsCluster.push_back(tdsLayer) ;
    }
}

}
