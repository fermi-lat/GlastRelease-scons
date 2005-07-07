
#include <RootConvert/Recon/CalParamsConvert.h>
#include <RootConvert/Utilities/Toolkit.h>

namespace RootPersistence {

void convert( const Event::CalParams & tdsParams, CalParams & rootParams ) {
    rootParams.init
      ( tdsParams.getEnergy(),tdsParams.getEnergyErr(),
        convert(tdsParams.getCentroid()),
        tdsParams.getxPosxPos(),tdsParams.getxPosyPos(),tdsParams.getxPoszPos(),
        tdsParams.getyPosyPos(),tdsParams.getyPoszPos(),tdsParams.getzPoszPos(),
        convert(tdsParams.getAxis()),
        tdsParams.getxDirxDir(),tdsParams.getxDiryDir(),tdsParams.getxDirzDir(),
        tdsParams.getyDiryDir(),tdsParams.getyDirzDir(),tdsParams.getzDirzDir() ) ;
 }
 
void convert( const CalParams & rootParams, Event::CalParams & tdsParams )
{
    const TVector3 & rootCentroid = rootParams.getCentroid() ;
    const TVector3 & rootAxis = rootParams.getAxis() ;
    tdsParams = Event::CalParams
      ( rootParams.getEnergy(), rootParams.getEnergyErr(),
        rootCentroid.X(),rootCentroid.Y(),rootCentroid.Z(),
        rootParams.getxPosxPos(), rootParams.getxPosyPos(), rootParams.getxPoszPos(),
        rootParams.getyPosyPos(), rootParams.getyPoszPos(), rootParams.getzPoszPos(),
        rootAxis.X(),rootAxis.Y(),rootAxis.Z(),
        rootParams.getxDirxDir(), rootParams.getxDiryDir(), rootParams.getxDirzDir(),
        rootParams.getyDiryDir(), rootParams.getyDirzDir(), rootParams.getzDirzDir() ) ;
}

}
