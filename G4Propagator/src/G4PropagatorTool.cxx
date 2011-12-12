// File and Version Information:
// $Header$
//
// Description: Service for particle transport management
//
// Author(s):
//      T.Usher

#include "G4PropagatorTool.h"

//static ToolFactory<G4PropagatorTool> g4prop_factory;
//const IToolFactory& G4PropagatorToolFactory = g4prop_factory;
DECLARE_TOOL_FACTORY(G4PropagatorTool);

IG4GeometrySvc* G4PropagatorTool::geometrySvc    = 0;
IPropagator*    G4PropagatorTool::propagatorTool = 0;

G4PropagatorTool::G4PropagatorTool(const std::string& type, const std::string& name, const IInterface* parent) :
  AlgTool(type, name, parent)
{
  // Purpose and Method:  Gaudi initialization routine. 
  // Inputs:  None
  // Outputs:  Gaudi StatusCode
  // Dependencies: None
  // Restrictions None 

  //Declare additional interface
  declareInterface<IPropagatorTool>(this);
}

StatusCode G4PropagatorTool::initialize()
{
  MsgStream log(msgSvc(), name());

  log << MSG::INFO << "G4 Propagator Tool initializing now" << endreq;

  IService* iService = 0;
  if( serviceLocator()->getService( "G4GeometrySvc", iService, true).isFailure() ) 
  {
    log << MSG::ERROR << "Could not find the G4 Geometry Service!" << endreq;
    return StatusCode::FAILURE;
  }

  IG4GeometrySvc* gsv = dynamic_cast<IG4GeometrySvc*>(iService);

  geometrySvc = gsv;

  StatusCode sc = toolSvc()->retrieveTool("G4PropagationTool", propagatorTool);
  if (sc.isSuccess()) {log << MSG::INFO << "G4 Propagator Tool ready" << endreq;}
  else  { log << MSG::ERROR << "Couldn't retrieve G4PropagationTool " << endreq; }

  return sc;
}

G4PropagatorTool::~G4PropagatorTool()
{
}
