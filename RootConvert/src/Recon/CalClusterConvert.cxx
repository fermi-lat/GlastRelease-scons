
#include <RootConvert/Recon/CalClusterConvert.h>
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
      
    const Event::CalParams & tdsParams = tdsCluster.getCalParams() ;
    CalParams rootParams
      ( tdsParams.getEnergy(),tdsParams.getEnergyErr(),
        convert(tdsParams.getCentroid()),
        tdsParams.getxPosxPos(),tdsParams.getxPosyPos(),tdsParams.getxPoszPos(),
        tdsParams.getyPosyPos(),tdsParams.getyPoszPos(),tdsParams.getzPoszPos(),
        convert(tdsParams.getAxis()),
        tdsParams.getxDirxDir(),tdsParams.getxDiryDir(),tdsParams.getxDirzDir(),
        tdsParams.getyDiryDir(),tdsParams.getyDirzDir(),tdsParams.getzDirzDir() ) ;
      
    rootCluster.init
      ( rootLayers, rootParams,
        tdsCluster.getRmsLong(),
        tdsCluster.getRmsLongAsym(),
        tdsCluster.getRmsTrans(),
        (Int_t)tdsCluster.getNumTruncXtals(),
        tdsCluster.getStatusBits() ) ;
 }
 
void convert( const CalCluster & rootCluster, Event::CalCluster & tdsCluster )
{
    const CalParams & rootParams = rootCluster.getParams() ;
    const TVector3 & rootCentroid = rootParams.getCentroid() ;
    const TVector3 & rootAxis = rootParams.getAxis() ;
    Event::CalParams tdsParams
      ( rootParams.getEnergy(), rootParams.getEnergyErr(),
        rootCentroid.X(),rootCentroid.Y(),rootCentroid.Z(),
        rootParams.getxPosxPos(), rootParams.getxPosyPos(), rootParams.getxPoszPos(),
        rootParams.getyPosyPos(), rootParams.getyPoszPos(), rootParams.getzPoszPos(),
        rootAxis.X(),rootAxis.Y(),rootAxis.Z(),
        rootParams.getxDirxDir(), rootParams.getxDiryDir(), rootParams.getxDirzDir(),
        rootParams.getyDiryDir(), rootParams.getyDirzDir(), rootParams.getzDirzDir() ) ;
    
    tdsCluster.initialize
      ( tdsParams,
        rootCluster.getRmsLong(),
        rootCluster.getRmsTrans(),
        rootCluster.getRmsLongAsym(),
        rootCluster.getNumTruncXtals() ) ;
    
    UInt_t rootStatusBits = rootCluster.getStatusBits() ;
    Event::CalCluster::StatusBits tdsNotStatusBits
        = (Event::CalCluster::StatusBits)(~rootStatusBits) ;
    tdsCluster.clearStatusBit(tdsNotStatusBits) ;
    tdsCluster.setStatusBit(rootStatusBits) ;

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
