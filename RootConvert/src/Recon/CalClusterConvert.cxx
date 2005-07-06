
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
//        Double_t energyRoot = tdsLayer->getEnergy() ;
//        const Point & positionTds = tdsLayer->getPosition() ;
//        const Vector & rmsSpreadTds = tdsLayer->getRmsSpread() ;
//        TVector3 positionRoot(positionTds.x(),positionTds.y(),positionTds.z()) ;
//        TVector3 rmsSpreadRoot(rmsSpreadTds.x(),rmsSpreadTds.y(),rmsSpreadTds.z()) ;
        CalClusterLayerData layerRoot
          ( tdsLayer->getEnergy(),
            convert(tdsLayer->getPosition()),
            convert(tdsLayer->getRmsSpread())) ;
        rootLayers.push_back(layerRoot) ;
    }
      
    const Event::CalParams & tdsParams = tdsCluster.getCalParams() ;
//    const Point & centroidTds = tdsParams.getCentroid() ;
//    const Vector & axisTds = tdsParams.getAxis() ;
//    TVector3 centroidRoot(centroidTds.x(),centroidTds.y(),centroidTds.z()) ;
//    TVector3 axisRoot(axisTds.x(),axisTds.y(),axisTds.z()) ;
    CalParams rootParams
      ( tdsParams.getEnergy(),tdsParams.getEnergyErr(),
        convert(tdsParams.getCentroid()),
        tdsParams.getxPosxPos(),tdsParams.getxPosyPos(),tdsParams.getxPoszPos(),
        tdsParams.getyPosyPos(),tdsParams.getyPoszPos(),tdsParams.getzPoszPos(),
        convert(tdsParams.getAxis()),
        tdsParams.getxDirxDir(),tdsParams.getxDiryDir(),tdsParams.getxDirzDir(),
        tdsParams.getyDiryDir(),tdsParams.getyDirzDir(),tdsParams.getzDirzDir() ) ;
      
//    Double_t rootRmsLong = tdsCluster.getRmsLong() ;
//    Double_t rootRmsLongAsym = tdsCluster.getRmsLongAsym() ;
//    Double_t rootRmsTrans = tdsCluster.getRmsTrans() ;
//    Int_t rootNumTruncXtals = (Int_t)tdsCluster.getNumTruncXtals() ;
//    UInt_t rootStatusBits = tdsCluster.getStatusBits() ;
//    
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
//        const TVector3 & rootPosition = rootLayer.getPosition() ;
//        const TVector3 & rootRmsSpread = rootLayer.getRmsSpread() ;
//        Point tdsPosition(rootPosition.X(),rootPosition.Y(),rootPosition.Z()) ;
//        Vector tdsRmsSpread(rootRmsSpread.X(),rootRmsSpread.Y(),rootRmsSpread.Z()) ;
        Event::CalClusterLayerData tdsLayer
          ( rootLayer.getEnergy(),
            convertPoint(rootLayer.getPosition()),
            convertVector(rootLayer.getRmsSpread()) ) ;
        tdsCluster.push_back(tdsLayer) ;
    }
}

void random( CalCluster & cl, UInt_t /* i */, Float_t /* randNum */ ) {
 
    CalParams p(1.0, 2.0, TVector3(3.0, 4.0, 5.0), 6.0, 7.0, 8.0,
        9.0, 10.0, 11.0, TVector3(12., 13., 14.), 15., 16.0,
        17.0, 18.0, 19.0, 20.0);
    std::vector<CalClusterLayerData> clusLayerData ;
    UInt_t iclusLayer;
    for ( iclusLayer = 0; iclusLayer<ROOT_NUMCALLAYERS ; ++iclusLayer ) {
        CalClusterLayerData layer
           (1.0,TVector3(2.0,3.0,4.0),
            TVector3(5.0, 6.0, 7.0)) ;
        clusLayerData.push_back(layer) ;
    }
    cl.init(clusLayerData,p,1.0,2.0,3.0,4,5) ;
}

bool compare( const CalParams & cp1, const CalParams & cp2 ) {
    
    bool result = true ;
    result = result && compare(cp1.getEnergy(),cp2.getEnergy(),"Params Energy") ;
    result = result && compare(cp1.getEnergyErr(),cp2.getEnergyErr(),"Params EnergyErr") ;
    result = result && compare(cp1.getCentroid(),cp2.getCentroid(),"Params Centroid") ;
    result = result && compare(cp1.getxPosxPos(),cp2.getxPosxPos(),"Params xPosxPos") ;
    result = result && compare(cp1.getxPosyPos(),cp2.getxPosyPos(),"Params xPosyPos") ;
    result = result && compare(cp1.getxPoszPos(),cp2.getxPoszPos(),"Params xPoszPos") ;
    result = result && compare(cp1.getyPosyPos(),cp2.getyPosyPos(),"Params yPosyPos") ;
    result = result && compare(cp1.getyPoszPos(),cp2.getyPoszPos(),"Params yPoszPos") ;
    result = result && compare(cp1.getzPoszPos(),cp2.getzPoszPos(),"Params zPoszPos") ;
    result = result && compare(cp1.getAxis(),cp2.getAxis(),"Params Axis") ;
    result = result && compare(cp1.getxDirxDir(),cp2.getxDirxDir(),"Params xDirxDir") ;
    result = result && compare(cp1.getxDiryDir(),cp2.getxDiryDir(),"Params xDiryDir") ;
    result = result && compare(cp1.getxDirzDir(),cp2.getxDirzDir(),"Params xDirzDir") ;
    result = result && compare(cp1.getyDiryDir(),cp2.getyDiryDir(),"Params yDiryDir") ;
    result = result && compare(cp1.getyDirzDir(),cp2.getyDirzDir(),"Params yDirzDir") ;
    result = result && compare(cp1.getzDirzDir(),cp2.getzDirzDir(),"Params zDirzDir") ;
    return result ;
        
}

bool compare( const CalClusterLayerData & cl1, const CalClusterLayerData & cl2 ) {

    bool result = true ;
    result = result && compare(cl1.getEnergy(),cl2.getEnergy(),"Layer Energy") ;
    result = result && compare(cl1.getPosition(),cl2.getPosition(),"Layer Position") ;
    result = result && compare(cl1.getRmsSpread(),cl2.getRmsSpread(),"Layer RmsSpread") ;
    return result ;

}

bool compare( const CalCluster & c1, const CalCluster & c2 ) {
    bool result = true ;
    result = result && compare(c1.getParams(),c2.getParams()) ;
    result = result && compare(ROOT_NUMCALLAYERS,NUMCALLAYERS,"NUMCALLAYERS") ;
    int i ;
    for ( i=0 ; i<ROOT_NUMCALLAYERS ; ++i ) {
        result = result && compare(c1.getLayer(i),c2.getLayer(i)) ;
    }
    result = result && compare(c1.getRmsLong(),c2.getRmsLong(),"RmsLong") ;
    result = result && compare(c1.getRmsLongAsym(),c2.getRmsLongAsym(),"RmsLongAsym") ;
    result = result && compare(c1.getRmsTrans(),c2.getRmsTrans(),"RmsTrans") ;
    result = result && compare(c1.getNumTruncXtals(),c2.getNumTruncXtals(),"NumTruncXtals") ;
    result = result && compare(c1.getStatusBits(),c2.getStatusBits(),"StatusBits") ;
    return result ;
}

}
