
#include <RootConvert/Recon/CalXtalsParamsConvert.h>
#include <RootConvert/Utilities/Toolkit.h>

namespace RootPersistence {

void convert( const Event::CalXtalsParams & tdsXtalsParams, CalXtalsParams & rootXtalsParams ) {
    rootXtalsParams.init
      ( tdsXtalsParams.getNumXtals(),tdsXtalsParams.getNumTruncXtals(),
	tdsXtalsParams.getNumSaturatedXtals(),
	tdsXtalsParams.getXtalRawEneSum(),tdsXtalsParams.getXtalCorrEneSum(),
	tdsXtalsParams.getXtalEneMax(),
	tdsXtalsParams.getXtalEneRms(),tdsXtalsParams.getXtalEneSkewness(),
	convert(tdsXtalsParams.getCentroid()) ) ;
 }
 
void convert( const CalXtalsParams & rootXtalsParams, Event::CalXtalsParams & tdsXtalsParams )
{
    tdsXtalsParams = Event::CalXtalsParams
      ( rootXtalsParams.getNumXtals(),rootXtalsParams.getNumTruncXtals(),
	rootXtalsParams.getNumSaturatedXtals(),
	rootXtalsParams.getXtalRawEneSum(),rootXtalsParams.getXtalCorrEneSum(),
	rootXtalsParams.getXtalEneMax(),
	rootXtalsParams.getXtalEneRms(),rootXtalsParams.getXtalEneSkewness(),
	convertPoint(rootXtalsParams.getCentroid()) ) ;
}

}
