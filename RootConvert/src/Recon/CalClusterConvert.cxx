#include <RootConvert/Recon/CalClusterConvert.h>
#include <Event/Recon/CalRecon/CalCluster.h>   
#include <reconRootData/CalCluster.h>

namespace RootPersistence {

void convert( const Event::CalCluster & tdsCluster, CalCluster & rootCluster )
 {
 /*
      std::vector<CalClusterLayerData> rootLayers ;
      std::vector<Event::CalClusterLayerData>::const_iterator tdsLayer ;
      for ( tdsLayer=tdsCluster.begin() ;
            tdsLayer!=tdsCluster.end() ;
            ++tdsLayer )
       {
        Double_t rootEnergy = tdsLayer->getEnergy() ;
        const Point & tdsPosition = tdsLayer->getPosition() ;
        const Vector & tdsRmsSpread = tdsLayer->getRmsSpread() ;
        TVector3 rootPosition(tdsPosition.x(),tdsPosition.y(),tdsPosition.z()) ;
        TVector3 rootRmsSpread(tdsRmsSpread.x(),tdsRmsSpread.y(),tdsRmsSpread.z()) ;
        CalClusterLayerData rootLayer(rootEnergy,rootPosition,rootRmsSpread) ;
        rootLayers.push_back(rootLayer) ;
       }
      
      const Event::CalParams & tdsParams = tdsCluster.getCalParams() ;
      const Point & tdsCentroid = tdsParams.getCentroid() ;
      const Vector & tdsAxis = tdsParams.getAxis() ;
      TVector3 rootCentroid(tdsCentroid.x(),tdsCentroid.y(),tdsCentroid.z()) ;
      TVector3 rootAxis(tdsAxis.x(),tdsAxis.y(),tdsAxis.z()) ;
      CalParams rootParams
       ( tdsParams.getEnergy(),tdsParams.getEnergyErr(),
         rootCentroid,
         tdsParams.getxPosxPos(),tdsParams.getxPosyPos(),tdsParams.getxPoszPos(),
         tdsParams.getyPosyPos(),tdsParams.getyPoszPos(),tdsParams.getzPoszPos(),
         rootAxis,
         tdsParams.getxDirxDir(),tdsParams.getxDiryDir(),tdsParams.getxDirzDir(),
         tdsParams.getyDiryDir(),tdsParams.getyDirzDir(),tdsParams.getzDirzDir() ) ;
      
      Double_t rootRmsLong = tdsCluster.getRmsLong() ;
      Double_t rootRmsLongAsym = tdsCluster.getRmsLongAsym() ;
      Double_t rootRmsTrans = tdsCluster.getRmsTrans() ;
      Int_t rootNumTruncXtals = (Int_t)tdsCluster.getNumTruncXtals() ;
      UInt_t rootStatusBits = tdsCluster.getStatusBits() ;

      rootCluster.init(rootLayers,rootParams,
        rootRmsLong,rootRmsLongAsym,rootRmsTrans,
        rootNumTruncXtals,rootStatusBits) ;
*/
    std::vector<CalClusterLayerData> rootLayers ;
    std::vector<Event::CalClusterLayerData>::const_iterator tdsLayer ;
    for ( tdsLayer=tdsCluster.begin() ;
          tdsLayer!=tdsCluster.end() ;
          ++tdsLayer )
    {
        Double_t energyRoot = tdsLayer->getEnergy() ;
        const Point & positionTds = tdsLayer->getPosition() ;
        const Vector & rmsSpreadTds = tdsLayer->getRmsSpread() ;
        TVector3 positionRoot(positionTds.x(),positionTds.y(),positionTds.z()) ;
        TVector3 rmsSpreadRoot(rmsSpreadTds.x(),rmsSpreadTds.y(),rmsSpreadTds.z()) ;
        CalClusterLayerData layerRoot(energyRoot,positionRoot,rmsSpreadRoot) ;
        rootLayers.push_back(layerRoot) ;
    }
      
    const Event::CalParams & tdsParams = tdsCluster.getCalParams() ;
    const Point & centroidTds = tdsParams.getCentroid() ;
    const Vector & axisTds = tdsParams.getAxis() ;
    TVector3 centroidRoot(centroidTds.x(),centroidTds.y(),centroidTds.z()) ;
    TVector3 axisRoot(axisTds.x(),axisTds.y(),axisTds.z()) ;
    CalParams rootParams
      ( tdsParams.getEnergy(),tdsParams.getEnergyErr(),
        centroidRoot,
        tdsParams.getxPosxPos(),tdsParams.getxPosyPos(),tdsParams.getxPoszPos(),
        tdsParams.getyPosyPos(),tdsParams.getyPoszPos(),tdsParams.getzPoszPos(),
        axisRoot,
        tdsParams.getxDirxDir(),tdsParams.getxDiryDir(),tdsParams.getxDirzDir(),
        tdsParams.getyDiryDir(),tdsParams.getyDirzDir(),tdsParams.getzDirzDir() ) ;
      
    Double_t rootRmsLong = tdsCluster.getRmsLong() ;
    Double_t rootRmsLongAsym = tdsCluster.getRmsLongAsym() ;
    Double_t rootRmsTrans = tdsCluster.getRmsTrans() ;
    Int_t rootNumTruncXtals = (Int_t)tdsCluster.getNumTruncXtals() ;
    UInt_t rootStatusBits = tdsCluster.getStatusBits() ;
    
    rootCluster.init(rootLayers, rootParams, rootRmsLong,rootRmsLongAsym,rootRmsTrans,
        rootNumTruncXtals,rootStatusBits);
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
    
    tdsCluster.initialize( tdsParams,
                           rootCluster.getRmsLong(),
                           rootCluster.getRmsLongAsym(),
                           rootCluster.getRmsTrans(),
                           rootCluster.getNumTruncXtals() ) ;
    
    UInt_t rootStatusBits = rootCluster.getStatusBits() ;
    Event::CalCluster::StatusBits tdsNotStatusBits
        = (Event::CalCluster::StatusBits)(~rootStatusBits) ;
    tdsCluster.clearStatusBit(tdsNotStatusBits) ;
    tdsCluster.setStatusBit(rootStatusBits) ;

    int i ;
    for ( i=0 ; i<NUMCALLAYERS ; ++i )
    {
        const CalClusterLayerData & rootLayer = rootCluster.getLayer(i) ;
        const TVector3 & rootPosition = rootLayer.getPosition() ;
        const TVector3 & rootRmsSpread = rootLayer.getRmsSpread() ;
        Point tdsPosition(rootPosition.X(),rootPosition.Y(),rootPosition.Z()) ;
        Vector tdsRmsSpread(rootRmsSpread.X(),rootRmsSpread.Y(),rootRmsSpread.Z()) ;
        Event::CalClusterLayerData tdsLayer
            ( rootLayer.getEnergy(),tdsPosition,tdsRmsSpread ) ;
        tdsCluster.push_back(tdsLayer) ;
    }
}

}
