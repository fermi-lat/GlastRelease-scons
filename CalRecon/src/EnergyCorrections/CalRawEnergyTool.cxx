/** @file CalRawEnergyTool.cxx
@brief implementation of the class CalRawEnergyTool

$Header$

*/

#include <CalRecon/ICalEnergyCorr.h>
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/GaudiException.h" 

#include "Event/TopLevel/EventModel.h"
#include "Event/Recon/CalRecon/CalCluster.h"
#include "CLHEP/Matrix/Vector.h"

/**   
* @class CalRawEnergyTool
* @author CalRecon Rewrite Group
*
* This sets the "raw" energy for an event (after clustering)
*
* $Header$
*/

class CalRawEnergyTool : public AlgTool, virtual public ICalEnergyCorr 
{
public:

    //! Constructor
    CalRawEnergyTool( const std::string& type, const std::string& name, const IInterface* parent);
    //! destructor
    virtual ~CalRawEnergyTool() {}; 

    StatusCode initialize();

    // worker function to get the corrected energy      
    Event::CalCorToolResult* doEnergyCorr(Event::CalClusterCol*, Event::TkrVertex* );

private:

    /// Pointer to the Gaudi data provider service
    IDataProviderSvc* m_dataSvc;
};

#include "GaudiKernel/DeclareFactoryEntries.h"
DECLARE_TOOL_FACTORY(CalRawEnergyTool) ;

CalRawEnergyTool::CalRawEnergyTool( const std::string & type, 
                                    const std::string & name, 
                                    const IInterface * parent )
                                    : AlgTool(type,name,parent)
{ 
    declareInterface<ICalEnergyCorr>(this) ; 
}

StatusCode CalRawEnergyTool::initialize()
{
    // This function does following initialization actions:
    //    - Sets up a pointer to the data service as a convenience for the tool
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;
    log << MSG::INFO << "Initializing CalRawEnergyTool" <<endreq;

    //Locate and store a pointer to the data service which allows access to the TDS
    if ((sc = service("EventDataSvc", m_dataSvc)).isFailure())
    {
        throw GaudiException("Service [EventDataSvc] not found", name(), sc);
    }

    return sc;
}


Event::CalCorToolResult* CalRawEnergyTool::doEnergyCorr(Event::CalClusterCol* calClusters, Event::TkrVertex* )
{
    //Purpose and method:
    //
    //   This function determines the overal "raw" energy for an event
    //   Note: This sums the energy contribution from all clusters and
    //         does a weighted average for centroid position and axis
    // 
    // TDS input:  CalClusters
    // TDS output: CalEventEnergy
    
    MsgStream lm(msgSvc(), name());

    // Set up to loop over all clusters to get total raw energy
    double     rawEnergy   = 0.;
    double     rawEneError = 0.;
    CLHEP::HepVector  posSum(3);
    CLHEP::HepMatrix  posWghtSum(3,3,0.);
    CLHEP::HepVector  axisSum(3);
    CLHEP::HepMatrix  axisWghtSum(3,3,0.);

    // Do the loop and accumulate information
    for(Event::CalClusterCol::iterator clusIter = calClusters->begin(); clusIter != calClusters->end(); clusIter++)
    {
        // Note that if we are doing "clustering" (ie multiple clusters) then the first cluster is the "uber" cluster
        // So... right now, break out if after the first cluster
        if (clusIter != calClusters->begin()) break;

        Event::CalCluster*       cluster = *clusIter;
        const Event::CalParams&  params  = cluster->getCalParams();

        CLHEP::HepVector centroid(3);
        centroid[0] = params.getCentroid().x();
        centroid[1] = params.getCentroid().y();
        centroid[2] = params.getCentroid().z();

        CLHEP::HepVector axis(3);
        axis[0] = params.getAxis().x();
        axis[1] = params.getAxis().y();
        axis[2] = params.getAxis().z();

        rawEnergy   += params.getEnergy();
        rawEneError += params.getEnergyErr() * params.getEnergyErr();

        CLHEP::HepMatrix posCovInv = params.getCentroidErrs();
        int       matInvErr = 0;
        posCovInv.invert(matInvErr);
        posWghtSum += posCovInv;
        posSum     += posCovInv * centroid;


        CLHEP::HepMatrix axisCovInv = params.getAxisErrs();
        axisCovInv.invert(matInvErr);
        axisWghtSum += axisCovInv;
        axisSum     += axisCovInv * axis;
    }

    // Get new errors and weighted average centroid 
    int matInvErr = 0;
    posWghtSum.invert(matInvErr);
    posSum = posWghtSum * posSum;

    // Get new errors and weighted average axis
    axisWghtSum.invert(matInvErr);
    axisSum = axisWghtSum * axisSum;

    // New estimate of energy error
    rawEneError = sqrt(rawEneError);

    // Create a CalParams object to contain the results
    Point  centroid(posSum[0], posSum[1], posSum[2]);
    Vector axis(axisSum[0], axisSum[1], axisSum[2]);
    Event::CalParams params(rawEnergy, rawEneError, centroid, posWghtSum, axis, axisWghtSum);

    // Create a CalCorToolResult object to hold the information
    Event::CalCorToolResult* corResult = new Event::CalCorToolResult();

//    const std::string& nm = name();
//    const std::string& ty = type();

    corResult->setStatusBit(Event::CalCorToolResult::VALIDPARAMS);
    corResult->setCorrectionName(type());
    corResult->setParams(params);
    corResult->setChiSquare(1.);

    return corResult;
}

