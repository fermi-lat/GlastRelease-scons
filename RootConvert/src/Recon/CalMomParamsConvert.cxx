
#include <RootConvert/Recon/CalMomParamsConvert.h>
#include <RootConvert/Utilities/Toolkit.h>

namespace RootPersistence {

void convert( const Event::CalMomParams & tdsMomParams, CalMomParams & rootMomParams ) {
    rootMomParams.init
      ( tdsMomParams.getEnergy(),tdsMomParams.getEnergyErr(),
        convert(tdsMomParams.getCentroid()),
        tdsMomParams.getxPosxPos(),tdsMomParams.getxPosyPos(),tdsMomParams.getxPoszPos(),
        tdsMomParams.getyPosyPos(),tdsMomParams.getyPoszPos(),tdsMomParams.getzPoszPos(),
        convert(tdsMomParams.getAxis()),
        tdsMomParams.getxDirxDir(),tdsMomParams.getxDiryDir(),tdsMomParams.getxDirzDir(),
        tdsMomParams.getyDiryDir(),tdsMomParams.getyDirzDir(),tdsMomParams.getzDirzDir(),
	tdsMomParams.getNumIterations(),tdsMomParams.getNumCoreXtals(),tdsMomParams.getNumXtals(),
	tdsMomParams.getTransRms(),tdsMomParams.getLongRms(),tdsMomParams.getLongRmsAsym(),
	tdsMomParams.getLongSkewness(),tdsMomParams.getCoreEnergyFrac(),tdsMomParams.getdEdxAverage(),
	tdsMomParams.getdEdxSpread() ) ;
 }
 
void convert( const CalMomParams & rootMomParams, Event::CalMomParams & tdsMomParams )
{
    const TVector3 & rootCentroid = rootMomParams.getCentroid() ;
    const TVector3 & rootAxis = rootMomParams.getAxis() ;
    tdsMomParams = Event::CalMomParams
      ( rootMomParams.getEnergy(), rootMomParams.getEnergyErr(),
        rootCentroid.X(),rootCentroid.Y(),rootCentroid.Z(),
        rootMomParams.getxPosxPos(), rootMomParams.getxPosyPos(), rootMomParams.getxPoszPos(),
        rootMomParams.getyPosyPos(), rootMomParams.getyPoszPos(), rootMomParams.getzPoszPos(),
        rootAxis.X(),rootAxis.Y(),rootAxis.Z(),
        rootMomParams.getxDirxDir(), rootMomParams.getxDiryDir(), rootMomParams.getxDirzDir(),
        rootMomParams.getyDiryDir(), rootMomParams.getyDirzDir(), rootMomParams.getzDirzDir(),
	rootMomParams.getNumIterations(),rootMomParams.getNumCoreXtals(),rootMomParams.getNumXtals(),
	rootMomParams.getTransRms(),rootMomParams.getLongRms(),rootMomParams.getLongRmsAsym(),
	rootMomParams.getLongSkewness(),rootMomParams.getCoreEnergyFrac(),rootMomParams.getdEdxAverage(),
	rootMomParams.getdEdxSpread() ) ;
}

}
