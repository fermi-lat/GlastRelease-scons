
#include "CalValsCorrTool.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/ToolFactory.h"
#include "AnalysisNtuple/IValsTool.h"

static const ToolFactory<CalValsCorrTool>  s_factory;
const IToolFactory& CalValsCorrToolFactory = s_factory;


CalValsCorrTool::CalValsCorrTool( const std::string& type, 
                                     const std::string& name, 
                                     const IInterface* parent)
                                     : AlgTool(type,name,parent)
{
	// declare base interface for all consecutive concrete classes
	declareInterface<IEnergyCorr>(this);
    // Declare the properties that may be set in the job options file
    declareProperty ("calValsToolName", m_calValsToolName="CalValsTool");
	
}



StatusCode CalValsCorrTool::initialize()

// This function does following initialization actions:
//    - extracts geometry constants from xml file using GlastDetSvc

{
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;
	log << MSG::INFO << "Initializing CalValsCorrTool" <<endreq;

    // it appears that retrieveTool cannot store a EnergyCorr* pointer. Will need
    // to find IEnergyCorr* and downcast it to EnergyCorr*

    sc = toolSvc()->retrieveTool(m_calValsToolName,m_calValsTool);
    if (sc.isFailure() ) {
        log << MSG::ERROR << "  Unable to create " << m_calValsToolName << endreq;
        return sc;
    }

    return sc;
}


StatusCode CalValsCorrTool::doEnergyCorr(double eTotal, Event::CalCluster* cluster)

//Purpose and method:
//
//   This function performs the calorimeter cluster reconstruction.
//   The main actions are:
//      - calculate energy sum
//                  energy per layer
//                  average position per layer
//                  quadratic spread per layer
//      - fit the particle direction using Fit_Direction() function
//      - calculate particle energy by profile fitting method
//          using Profile() function
//      - calculate particle energy by last layer correction method
//          using Leak() function
//      - store all calculated quantities in CalCluster object
// 
// TDS input: CalXtalRecCol
// TDS output: CalClustersCol


{
    MsgStream lm(msgSvc(), name());
    StatusCode sc = m_calValsTool->calculate();
    
    double correctedEnergy;
    sc = m_calValsTool->getVal("CAL_Energy_Corr",correctedEnergy);
    if (sc == StatusCode::SUCCESS) cluster->setEnergyCorrected(correctedEnergy);

    return sc;
}

StatusCode CalValsCorrTool::finalize()
{
	StatusCode sc = StatusCode::SUCCESS;

	return sc;
}

StatusCode CalValsCorrTool::execute()
{
	StatusCode sc = StatusCode::SUCCESS;

	return sc;
}

