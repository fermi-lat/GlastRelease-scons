// File and Version Information:
// $Header$
//
// Description: Service for particle transport management
//
// Author(s):
//      T.Usher

#include "G4PropagatorTool.h"

static ToolFactory<G4PropagatorTool> g4prop_factory;
const IToolFactory& G4PropagatorToolFactory = g4prop_factory;

G4VUserDetectorConstruction* G4PropagatorTool::UserDetector          = 0;
IG4GeometrySvc::IdMap*       G4PropagatorTool::VolIdentMap           = 0;
G4TransportationManager*     G4PropagatorTool::TransportationManager = 0;

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

  //Look up message service
  IService* iService = 0;
  StatusCode sc = serviceLocator()->getService("MessageSvc", iService, true);

  if (sc.isSuccess())
  {
	IMessageSvc* msgSvc = dynamic_cast<IMessageSvc*>(iService);

	MsgStream log(msgSvc, name);

	log << MSG::DEBUG << "Propagator Tool initializing now" << endreq;

	if( serviceLocator()->getService( "G4GeometrySvc", iService).isFailure() ) 
	{
	    log << MSG::ERROR << "Could not find the G4 Geometry Service!" << endreq;
        return;
	}

	IG4GeometrySvc* gsv = dynamic_cast<IG4GeometrySvc*>(iService);

	UserDetector          = gsv->getDetector();
	TransportationManager = gsv->getTransportationManager();
	VolIdentMap           = gsv->getIdMap();

	log << MSG::DEBUG << "Propagator Tool ready" << endreq;
  }

  return;
}

G4PropagatorTool::~G4PropagatorTool()
{
}
