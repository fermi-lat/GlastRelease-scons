
#include <RootConvert/Recon/CalClusterConvert.h>
#include <RootConvert/Recon/CalFitParamsConvert.h>
#include <RootConvert/Recon/CalMSTreeParamsConvert.h>
#include <RootConvert/Recon/CalParamsConvert.h>
#include <RootConvert/Recon/CalMomParamsConvert.h>
#include <RootConvert/Recon/CalClassParamsConvert.h>
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
      
    CalMSTreeParams rootTreeParams ;
    convert(tdsCluster.getMSTreeParams(),rootTreeParams) ;
    CalFitParams rootFitParams ;
    convert(tdsCluster.getFitParams(),rootFitParams) ;
    CalMomParams rootMomParams ;
    convert(tdsCluster.getMomParams(),rootMomParams) ;
    CalClassParams rootClassParams;
    convert(tdsCluster.getClassParams(),rootClassParams) ;

    rootCluster.init
      ( rootLayers,
        rootTreeParams,
	rootFitParams,
	rootMomParams,
	rootClassParams,
        (Int_t)tdsCluster.getNumSaturatedXtals(),
        (Int_t)tdsCluster.getNumTruncXtals(),
        tdsCluster.getStatusBits() ) ;
 }
 
void convert( const CalCluster & rootCluster, Event::CalCluster & tdsCluster )
{
    Event::CalMSTreeParams tdsTreeParams;
    convert(rootCluster.getMSTreeParams(),tdsTreeParams) ;    
    Event::CalFitParams tdsFitParams;
    convert(rootCluster.getFitParams(),tdsFitParams) ;    
    Event::CalMomParams tdsMomParams ;
    convert(rootCluster.getMomParams(),tdsMomParams) ;
    Event::CalClassParams tdsClassParams ;
    convert(rootCluster.getClassParams(),tdsClassParams) ;

    tdsCluster.initialize
      ( tdsTreeParams, 
        tdsFitParams,
        tdsMomParams,
	tdsClassParams,
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
