
#include <RootConvert/Recon/CalMSTreeParamsConvert.h>
#include <RootConvert/Utilities/Toolkit.h>

namespace RootPersistence {

void convert( const Event::CalMSTreeParams & tdsParams, CalMSTreeParams & rootParams )
{
    rootParams.init
      ( tdsParams.getTotalEnergy(),
        tdsParams.getMaxXtalEnergy(),tdsParams.getNumberOfEdges(),
	tdsParams.getMinEdgeLength(),tdsParams.getMaxEdgeLength(),
	tdsParams.getMeanEdgeLength(),tdsParams.getMeanEdgeLengthTrunc(),
	tdsParams.getRmsEdgeLength(),tdsParams.getRmsEdgeLengthTrunc() );
 }
 
void convert( const CalMSTreeParams & rootParams, Event::CalMSTreeParams & tdsParams )
{
    tdsParams = Event::CalMSTreeParams
      ( rootParams.getTotalEnergy(),
        rootParams.getMaxXtalEnergy(),rootParams.getNumberOfEdges(),
	rootParams.getMinEdgeLength(),rootParams.getMaxEdgeLength(),
	rootParams.getMeanEdgeLength(),rootParams.getMeanEdgeLengthTrunc(),
	rootParams.getRmsEdgeLength(),rootParams.getRmsEdgeLengthTrunc()) ;
}

}
