
/** @file EvtValsTool.cxx
@brief Calculates the "Event" analysis variables from the other ntuple variables
@author Bill Atwood, Leon Rochester

$Header$
*/

#include "ValBase.h"

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/IToolSvc.h"

#include "Event/TopLevel/EventModel.h"
#include "Event/TopLevel/Event.h"

#include "TkrUtil/ITkrGeometrySvc.h"

#include <algorithm>
/** @class EvtValsTool
@brief Calculates Event Values from the other ntuple variables
*/

class EvtValsTool : public ValBase
{
public:

    EvtValsTool( const std::string& type, 
        const std::string& name, 
        const IInterface* parent);

    virtual ~EvtValsTool() { }

    StatusCode initialize();

    StatusCode calculate();

private:

    //Global Items
    float EvtRun;
    float EvtEventId;
    double EvtElapsedTime;
    float EvtLiveTime;

    float EvtEnergyCorr;
    float EvtEnergyRaw;
    float EvtDeltaEoE;
    float EvtCalEdgeAngle;
    float EvtTkrEdgeAngle;
    float EvtLogEnergy;
    float EvtTkr1EFrac;
    float EvtVtxKin;
    float EvtVtxEAngle;
    float EvtTkrComptonRatio;
    float EvtETkrComptonRatio;
    float EvtPSFModel; 

    float EvtETkr1Chisq;
    float EvtETkr1FirstChisq;
    float EvtETkr1Qual;
    float EvtTkr1PSFMdRat;

    float EvtECalXtalRatio;
    float EvtECalXtalTrunc;
    float EvtECalTrackDoca;
    float EvtECalTrackSep;
	float EvtECalTransRms;
	float EvtECalLongRms;
	float EvtECalLRmsAsym;
	float EvtECalTrackAngle;
    float EvtEVtxAngle;
    float EvtEVtxDoca;

    IValsTool* m_pMcTool;
    IValsTool* m_pGltTool;
    IValsTool* m_pTkrHitTool;
    IValsTool* m_pTkrTool;
    IValsTool* m_pVtxTool;
    IValsTool* m_pCalTool;
    IValsTool* m_pAcdTool;

    ITkrGeometrySvc* m_tkrGeom;

	// These are copies of the new DC2 vars from IM.   Left here for reference
// EvtECalXtalRatio	continuous	CalXtalRatio/(2.99-1.19*LogEnergyOpt  + .122*LogEnergyOpt^2)/(.749-.355*Tkr1ZDir)
// EvtECalTrackDoca	continuous	CalTrackDoca/(272-140.5*min(EvtLogEnergy, 3.8) + 18.7*min(EvtLogEnergy, 3.8)^2)/(3.08+2.67*Tkr1ZDir)
// EvtECalXtalsTrunc	continuous	CalXtalsTrunc/(max(1.,(-85.4 + 65.2*LogEnergyOpt - 10.4*LogEnergyOpt^2)) + min(12*(LogEnergyOpt-3.3)^2*(max(3.3,LogEnergyOpt)-3.3)/(LogEnergyOpt-3.3),14.))/(1.2+.224*Tkr1ZDir)
// EvtECalTransRms	continuous	CalTransRms/(12.5+ 22.*exp(-(EvtLogEnergy-2.3)^2/.8))/(1.34+.55*Tkr1ZDir)
// EvtECalLongRms	continuous	CalLongRms/(ifelse(EvtLogEnergy < 3.8, (96.1 - 39*EvtLogEnergy + 5.6*EvtLogEnergy^2),  28.*(.847-.0134*EvtLogEnergy+.015*EvtLogEnergy^2)))/(1.06+.0867*Tkr1ZDir)
// EvtECalLRmsAsym	continuous	CalLRmsAsym/(.012+.06*exp(-(EvtLogEnergy-2.0)^2/1.3))/(1.01-.0718*Tkr1ZDir)
// EvtECalTrackAngle	continuous	CalTrackAngle/(2.3-1.11*min(EvtLogEnergy,4.0)+.138*min(EvtLogEnergy,4.0)^2)/(1.43+.612*Tkr1ZDir)

// EvtETkr1Chisq	continuous	Tkr1Chisq/(6.49 -22.1*Tkr1ZDir -22.8*Tkr1ZDir^2)
// EvtETkr1FirstChisq	continuous	Tkr1FirstChisq/(4.34-.708*EvtLogEnergy+.19*EvtLogEnergy^2)/(.751-1.74*Tkr1ZDir-1.83*Tkr1ZDir^2)
// EvtETkr1Qual	continuous	Tkr1Qual/(56.6 +6.78*Tkr1ZDir + 8.23*Tkr1ZDir^2)
// EvtETkrComptonRatio	continuous	EvtTkrComptonRatio / (-9.34+7.22*EvtLogEnergy - .672*EvtLogEnergy^2)/(-.509 - 3.65*Tkr1ZDir - 2.03*Tkr1ZDir^2)

// EvtVtxEAngle	continuous	VtxAngle*sqrt(EvtEnergySumOpt)/(4.24 -1.98*EvtLogEnergy + .269*EvtLogEnergy^2)/(1.95+2.36*Tkr1ZDir+1.3*Tkr1ZDir^2)


};

// Static factory for instantiation of algtool objects
static ToolFactory<EvtValsTool> s_factory;
const IToolFactory& EvtValsToolFactory = s_factory;

// Standard Constructor
EvtValsTool::EvtValsTool(const std::string& type, 
                         const std::string& name, 
                         const IInterface* parent)
                         : ValBase( type, name, parent )
{    
    // Declare additional interface
    declareInterface<IValsTool>(this);
}

StatusCode EvtValsTool::initialize()
{
    StatusCode sc = StatusCode::SUCCESS;
    StatusCode fail = StatusCode::FAILURE;

    MsgStream log(msgSvc(), name());

    if( ValBase::initialize().isFailure()) return StatusCode::FAILURE;

    // get the services

    // try using one tool from another

    if(service( "TkrGeometrySvc", m_tkrGeom, true ).isFailure()) {
        log << MSG::ERROR << "Could not find TkrGeometrySvc" << endreq;
        return fail;
    }

    IToolSvc* pToolSvc = 0; 
    sc = service("ToolSvc", pToolSvc, true);
    if (!sc.isSuccess ()){
        log << MSG::ERROR << "Can't find ToolSvc, will quit now" << endreq;
        return StatusCode::FAILURE;
    }

    m_pMcTool = 0;
    sc = pToolSvc->retrieveTool("McValsTool", m_pMcTool);
    if( sc.isFailure() ) {
        log << MSG::ERROR << "Unable to find tool: " "McValsTool" << endreq;
        return sc;
    }
    m_pGltTool = 0;
    sc = pToolSvc->retrieveTool("GltValsTool", m_pGltTool);
    if( sc.isFailure() ) {
        log << MSG::ERROR << "Unable to find tool: " "GltValsTool" << endreq;
        return sc;
    }
    m_pTkrHitTool = 0;
    sc = pToolSvc->retrieveTool("TkrHitValsTool", m_pTkrHitTool);
    if( sc.isFailure() ) {
        log << MSG::ERROR << "Unable to find tool: " "TkrHitValsTool" << endreq;
        return sc;
    }
    m_pTkrTool = 0;
    sc = pToolSvc->retrieveTool("TkrValsTool", m_pTkrTool);
    if( sc.isFailure() ) {
        log << MSG::ERROR << "Unable to find tool: " "TkrValsTool" << endreq;
        return sc;
    }
    m_pVtxTool = 0;
    sc = pToolSvc->retrieveTool("VtxValsTool", m_pVtxTool);
    if( sc.isFailure() ) {
        log << MSG::ERROR << "Unable to find tool: " "VtxValsTool" << endreq;
        return sc;
    }
    m_pCalTool = 0;
    sc = pToolSvc->retrieveTool("CalValsTool", m_pCalTool);
    if( sc.isFailure() ) {
        log << MSG::ERROR << "Unable to find tool: " "CalValsTool" << endreq;
        return sc;
    }
    m_pAcdTool = 0;
    sc = pToolSvc->retrieveTool("AcdValsTool", m_pAcdTool);
    if( sc.isFailure() ) {
        log << MSG::ERROR << "Unable to find tool: " "AcdValsTool" << endreq;
        return sc;
    }


    // load up the map
    addItem("EvtRun",           &EvtRun);
    addItem("EvtEventId",       &EvtEventId);
    addItem("EvtElapsedTime",   &EvtElapsedTime);
    addItem("EvtLiveTime",      &EvtLiveTime);

    addItem("EvtEnergyCorr",    &EvtEnergyCorr);
    addItem("EvtEnergyRaw",     &EvtEnergyRaw);
    addItem("EvtDeltaEoE",      &EvtDeltaEoE);
    addItem("EvtCalEdgeAngle",  &EvtCalEdgeAngle);
    addItem("EvtTkrEdgeAngle",  &EvtTkrEdgeAngle);
    addItem("EvtLogEnergy",     &EvtLogEnergy);
    addItem("EvtTkr1EFrac",     &EvtTkr1EFrac);
    addItem("EvtVtxKin",        &EvtVtxKin);
    addItem("EvtVtxEAngle",     &EvtVtxEAngle);
    addItem("EvtTkrComptonRatio", &EvtTkrComptonRatio);
	addItem("EvtETkrComptonRatio", &EvtETkrComptonRatio);
    addItem("EvtPSFModel",      &EvtPSFModel);

    addItem("EvtETkr1Chisq",    &EvtETkr1Chisq);
    addItem("EvtETkr1FirstChisq", &EvtETkr1FirstChisq);
    addItem("EvtETkr1Qual",     &EvtETkr1Qual);
    addItem("EvtTkr1PSFMdRat",  &EvtTkr1PSFMdRat);

    addItem("EvtECalXtalRatio", &EvtECalXtalRatio);
    addItem("EvtECalXtalTrunc", &EvtECalXtalTrunc);
    addItem("EvtECalTrackDoca", &EvtECalTrackDoca);
    addItem("EvtECalTrackSep",  &EvtECalTrackSep);
	addItem("EvtECalTransRms",  &EvtECalTransRms);
    addItem("EvtECalLongRms",   &EvtECalLongRms);
	addItem("EvtECalLRmsAsym",  &EvtECalLRmsAsym);
	addItem("EvtECalTrackAngle",&EvtECalTrackAngle);

    addItem("EvtEVtxAngle",    &EvtEVtxAngle);
    addItem("EvtEVtxDoca",      &EvtEVtxDoca);

    zeroVals();

    return sc;
}

StatusCode EvtValsTool::calculate()
{
    StatusCode sc = StatusCode::SUCCESS;

    MsgStream log(msgSvc(), name());

    // we may want to add TDS stuff to this method, but we haven't needed it yet.

    int firstCheck = m_check;
    int nextCheck = -1;

    // since we know what's happening here, we can plan a little
    // the idea is to call the first check of each tool with the called check value,
    // and the rest with the no-calc value
    // so be careful when adding calls or moving stuff around!!!!

    SmartDataPtr<Event::EventHeader> header(m_pEventSvc, EventModel::EventHeader);

    if(header) {
        EvtRun         = header->run();
        EvtEventId     = header->event();
        EvtElapsedTime = header->time();
        EvtLiveTime    = header->livetime();
    }

    double eCalSum, eTkr;
    if(    m_pCalTool->getVal("CalEnergyRaw", eCalSum, firstCheck).isSuccess()
        && m_pTkrTool->getVal("TkrEnergyCorr", eTkr, firstCheck).isSuccess()) {
        EvtEnergyRaw = eTkr + eCalSum;
    }

    double eTkrKalEne, eCalRLn, eTkrBest; //, eCal
    int CAL_Type;
    if (   m_pCalTool->getVal("CalCsIRLn", eCalRLn, firstCheck).isSuccess()
        && m_pTkrTool->getVal("TkrSumKalEne", eTkrKalEne, firstCheck).isSuccess()) 
    {
        eTkrBest = std::max(eTkr+eCalSum, eTkrKalEne);
        if(eCalSum < 350 || eCalRLn < 4) {
            if(eCalSum < 5 || eCalRLn < 4) {
                CAL_Type = 0;
            }
            else {CAL_Type = 1;}
        }
        else {CAL_Type = 2;} 
    }

 // EdgeAngle:   distance from tower edge / cos(theta) 
    double tkrEdge, calEdge, tkr1ZDir = -1.;
    if(m_pTkrTool->getVal("Tkr1ZDir",tkr1ZDir, nextCheck).isSuccess()) {
		if(tkr1ZDir == 0.) tkr1ZDir = -1.; 
        if (m_pTkrTool->getVal("TkrTwrEdge", tkrEdge, nextCheck).isSuccess()) {
            EvtTkrEdgeAngle = -tkrEdge/tkr1ZDir;
        }
        if (m_pCalTool->getVal("CalTwrEdge", calEdge, nextCheck).isSuccess()) {
            EvtCalEdgeAngle = -calEdge/tkr1ZDir;
        }
    }
	double tkr1ZDir2 = tkr1ZDir*tkr1ZDir;  
	// Final Energy Selection: (Combo approach) 
	//    Use Kalman energy when not enough in the Cal,
	//    use Last Layer Correction when avaialable (eCalEneLLCorr > 0) or else
	//    use integrated edge and leakage correction from CalValsTool
    double eCalSumCorr;
    if(    m_pCalTool->getVal("CalEnergyCorr", eCalSumCorr, nextCheck).isSuccess()
		//&& m_pCalTool->getVal("CalEnergyLLCorr", eCalEneLLCorr, nextCheck).isSuccess()
		) { // NOTE: THIS IS STILL LARGELY UNDESIDED WHAT TO DO HERE....!!!!!!!  
    //    if(CAL_Type == 0) EvtEnergyCorr = eTkrBest;
	//	else {
		//	if(eCalEneLLCorr > 0) EvtEnergyCorr = eCalEneLLCorr; //This is wrong! LL has comp.for TkrEne!
		//	else { 
				EvtEnergyCorr = (eTkr + eCalSumCorr);
		//	}
	//	}
    }
    
    double mcEnergy;
    if(m_pMcTool->getVal("McEnergy", mcEnergy, firstCheck).isSuccess()){
        EvtDeltaEoE = (EvtEnergyCorr - mcEnergy)/(mcEnergy);
    }
    // Model simple for PSF(68%) 
    EvtPSFModel = sqrt(pow((.061/pow((std::max(EvtEnergyCorr*1.,1.)/100),.8)),2) + (.001745*.001745));

	// Log(base 10) of measured energy - useful for parameterizing effects
    EvtLogEnergy = log10(std::min(std::max(EvtEnergyCorr,10.f),1000000.f));
    double logE = std::min(std::max(EvtLogEnergy,1.3f), 6.0f);
    double logE2 = logE*logE; 
    

	// Ratio of PSF model to event total PSF error. Note PhiErr = sin(theta)*d(phi)
    double tkr1ThetaErr, tkr1PhiErr;
    if (m_pTkrTool->getVal("Tkr1ThetaErr",tkr1ThetaErr, nextCheck).isSuccess() &&
        m_pTkrTool->getVal("Tkr1PhiErr",tkr1PhiErr, nextCheck).isSuccess()) {
        EvtTkr1PSFMdRat = sqrt(tkr1ThetaErr*tkr1ThetaErr + tkr1PhiErr*tkr1PhiErr)/ 
                          EvtPSFModel;
    }

	// Fraction of energy in Track 1
    double tkr1ConE;
    if (m_pTkrTool->getVal("Tkr1ConEne",tkr1ConE, nextCheck).isSuccess()) {
        if(EvtEnergyCorr>0.0) EvtTkr1EFrac = tkr1ConE/EvtEnergyCorr;
    }

	// Vtx kinematic variable:  angle * event energy / Track_1 energy fraction
    double vtxAngle;
    if (m_pVtxTool->getVal("VtxAngle", vtxAngle, firstCheck).isSuccess()) {
        if (tkr1ConE>0.0) EvtVtxKin = vtxAngle*EvtEnergyCorr/EvtTkr1EFrac;
    }

	// Vtx angle x event energy  ~ constant
    EvtVtxEAngle = vtxAngle*EvtEnergyCorr;
    
	// Hit counting around track compared to hits on track.  Compton scatters
	// have a ratio ~< 1 (as do MIPs etc.  - this is similar to the former
	// Surplus Hit Ratio
    double totHits, tkr1First;
    if (m_pTkrTool->getVal("TkrTotalHits", totHits, nextCheck).isSuccess()) {
        if (m_pTkrTool->getVal("Tkr1FirstLayer", tkr1First, nextCheck).isSuccess()){
            EvtTkrComptonRatio = totHits/(2.*(m_tkrGeom->numLayers()-tkr1First));
            EvtETkrComptonRatio = EvtTkrComptonRatio/(-9.34+7.22*logE - .672*logE2)
				                         /(-.509 - 3.65*tkr1ZDir - 2.03*tkr1ZDir2);
        }
    }

	// Energy compensated track 1 chisq
    double tkr1Chisq;
    if (m_pTkrTool->getVal("Tkr1Chisq", tkr1Chisq, nextCheck).isSuccess()) {
        EvtETkr1Chisq = tkr1Chisq/(6.49 -22.1*tkr1ZDir -22.8*tkr1ZDir2);
    }
    double tkr1_1stChisq;
    if (m_pTkrTool->getVal("Tkr1FirstChisq", tkr1_1stChisq, nextCheck).isSuccess()) {
        EvtETkr1FirstChisq = tkr1_1stChisq/(4.34-.708*logE+.19*logE2)
			                 /(.751-1.74*tkr1ZDir-1.83*tkr1ZDir2);
    }
    double tkr1Qual;
    if (m_pTkrTool->getVal("Tkr1Qual", tkr1Qual, nextCheck).isSuccess()) {
        EvtETkr1Qual = tkr1Qual/(56.6 +6.78*tkr1ZDir + 8.23*tkr1ZDir2);
    }

    double calXtalRatio;
    if(m_pCalTool->getVal("CalXtalRatio", calXtalRatio, nextCheck).isSuccess()) {
        EvtECalXtalRatio = calXtalRatio/(2.99-1.19*logE  + .122*logE2)/(.749-.355*tkr1ZDir);
    }

    double calXtalsTrunc;
    if(m_pCalTool->getVal("CalXtalsTrunc", calXtalsTrunc, nextCheck).isSuccess()) {
        double term1 = std::max(1., (-85.4 + 65.2*logE - 10.4*logE2));
        double logE33 = logE-3.3; 
        double term2 = (logE33 > 0.) ? std::min(14., 12.*logE33*logE33): 0.;
        EvtECalXtalTrunc = calXtalsTrunc/(term1 + term2)/(.935 - .382*tkr1ZDir - .343*tkr1ZDir2);
    }

    double calTrackDoca;
    if(m_pCalTool->getVal("CalTrackDoca", calTrackDoca, nextCheck).isSuccess()) {
		double logETrunc = std::min(EvtLogEnergy, 3.8f); 
        EvtECalTrackDoca = calTrackDoca/(272.-140.5*logETrunc + 18.7*logETrunc*logETrunc)
			                           /(3.08+2.67*tkr1ZDir);
    }

	double calTransRms;
    if(m_pCalTool->getVal("CalTransRms", calTransRms, nextCheck).isSuccess()) {
		double logEGaussSq = (EvtLogEnergy-2.3)*(EvtLogEnergy-2.3);
        EvtECalTransRms = calTransRms/(12.5+ 22.*exp(-logEGaussSq/.8))/(1.34+.55*tkr1ZDir);
    }

	double calLongRms;
    if(m_pCalTool->getVal("CalLongRms", calLongRms, nextCheck).isSuccess()) {
		if(logE < 3.8) {
		    EvtECalLongRms = calLongRms/(96.1 - 39*logE + 5.6*logE2);
		}
		else {
			EvtECalLongRms = calLongRms/(28.*(.847-.0134*logE+.015*logE2));
		}
		EvtECalLongRms /= (1.06+.0867*tkr1ZDir);
    }

	double calLRmsAsym;
    if(m_pCalTool->getVal("CalLRmsAsym", calLRmsAsym, nextCheck).isSuccess()) {
		double logEGaussSq = (logE-2.0)*(logE-2.0); 
		EvtECalLRmsAsym	 = calLRmsAsym/(.012+.06*exp(-logEGaussSq/1.3))/(1.01-.0718*tkr1ZDir);
    }

    double calTrackAngle;
    if(m_pCalTool->getVal("CalTrackAngle", calTrackAngle, nextCheck).isSuccess()) {
		double logETrunc = std::min(logE, 4.0); 
		EvtECalTrackAngle =	calTrackAngle/(2.3-1.11*logETrunc +.138*logETrunc*logETrunc)
			                /(1.43+.612*tkr1ZDir);
    }

    EvtEVtxAngle = vtxAngle*sqrt(EvtEnergyCorr)/(4.24 -1.98*logE + .269*logE2)
		                                       /(1.95+2.36*tkr1ZDir+1.3*tkr1ZDir2);

    double vtxDoca;
    if (m_pVtxTool->getVal("VtxDOCA", vtxDoca, firstCheck).isSuccess()) {
        EvtEVtxDoca = vtxDoca/(1.55 - .685*logE+ .0851*logE2) 
                                / (2.21 + 3.01*tkr1ZDir + 1.59*tkr1ZDir2);
    }

    return sc;
}
