
#include "CalValsCorrTool.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/ToolFactory.h"
#include "AnalysisNtuple/IValsTool.h"

static const ToolFactory<CalValsCorrTool>  s_factory;
const IToolFactory& CalValsCorrToolFactory = s_factory;


CalValsCorrTool::CalValsCorrTool( const std::string& type, 
                                     const std::string& name, 
                                     const IInterface* parent)
                                     : EnergyCorr(type,name,parent)
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


    sc = toolSvc()->retrieveTool(m_calValsToolName,m_calValsTool);
    if (sc.isFailure() ) {
        log << MSG::ERROR << "  Unable to create " << m_calValsToolName << endreq;
        return sc;
    }

    return sc;
}


StatusCode CalValsCorrTool::doEnergyCorr(double, Event::CalCluster* cluster)

//Purpose and method:
//
//   This function calls CalValsTool and extracts the crack/leakage corrected
//   energy, storing it in CalCluster.
// 
// TDS input: none
// TDS output: CalClusters


{
    MsgStream lm(msgSvc(), name());
    
    double correctedEnergy;
    StatusCode sc = m_calValsTool->getVal("CalEnergyCorr",correctedEnergy);
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

