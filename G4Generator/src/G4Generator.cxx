// $Header$

// Include files

// Geant4
#include "G4Generator.h"
#include "G4UImanager.hh"

#include "RunManager.h"
#include "PrimaryGeneratorAction.h"

// Gaudi

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/IAlgManager.h"
#include "GaudiKernel/ISvcLocator.h"

//flux
#include "FluxSvc/FluxSvc.h"
#include "FluxSvc/IFlux.h"

#include "CLHEP/Geometry/Point3D.h"
#include "CLHEP/Geometry/Vector3D.h"

static const AlgFactory<G4Generator>  Factory;
const IAlgFactory& G4GeneratorFactory = Factory;

G4Generator::G4Generator(const std::string& name, ISvcLocator* pSvcLocator) 
:Algorithm(name, pSvcLocator) 
{
// set defined properties
//    setProperty("file_name", m_file_name);
     declareProperty("source_name",  m_source_name="default");
     declareProperty("UIcommands", m_UIcommands);
     declareProperty("topVolume", m_topvol="");
     declareProperty("visitorMode", m_visitorMode="fastmc");

}
    
////////////////////////////////////////////////////////////////////////////
StatusCode G4Generator::initialize()
{
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "initialize" << endreq;

    // Use the Job options service to set the Algorithm's parameters
    setProperties();

    if ( service("FluxSvc", m_fluxSvc).isFailure() ){
        log << MSG::ERROR << "Couldn't find the FluxSvc!" << endreq;
        return StatusCode::FAILURE;
    }
    
    if ( m_fluxSvc->source(m_source_name, m_flux).isFailure() ){
        log << MSG::ERROR << "Couldn't find the source \"" 
            << m_source_name << "\"" << endreq;
        return StatusCode::FAILURE;
    }
    log << MSG::INFO << "Source: "<< m_flux->title() << endreq;

    // Set the geant4 classes needed for the simulation
    // The manager, with specified (maybe) top volume
    m_runManager = new RunManager(m_topvol, m_visitorMode);

    if( !m_UIcommands.value().empty() ) {
        G4UImanager* UI = G4UImanager::GetUIpointer();
        for( std::vector<std::string>::const_iterator k = m_UIcommands.value().begin(); 
            k!=m_UIcommands.value().end(); ++k){
            UI->ApplyCommand(*k);
            log << MSG::INFO << "UI command: " << (*k) << endreq;
        }
            UI->ApplyCommand("/event/verbose 3");
    }  
    // Initialize Geant4
    m_runManager->Initialize();

    return StatusCode::SUCCESS;
}
//------------------------------------------------------------------------------
StatusCode G4Generator::execute() 
{
    MsgStream   log( msgSvc(), name() );

    // following model from previous version, allow property "UIcommands" to generate
    // UI commands here. (but does not seem to work??)
    //
    if( !m_UIcommands.value().empty() ) {
        for( std::vector<std::string>::const_iterator k = m_UIcommands.value().begin(); 
            k!=m_UIcommands.value().end(); ++k){
            G4UImanager::GetUIpointer()->ApplyCommand(*k);
            log << MSG::INFO << "Apply UI command: \"" << (*k) << "\"" <<endreq;
        }
    }  

    //
    // have the flux service create parameters of an incoming particle, 
    // and define it as a MCParticle
    //
    m_flux->generate();

    // these are the particle properties
    std::string name(m_flux->particleName());
    HepVector3D dir(m_flux->launchDir());
    double ke= m_flux->energy() ;
    HepPoint3D p(m_flux->launchPoint());
    
    p = 10*p;
    ke = ke*1000;
    
    PrimaryGeneratorAction* primaryGenerator = 
      (PrimaryGeneratorAction*)m_runManager->GetUserPrimaryGeneratorAction();
    
    // Set the G4 primary generator
    // the position has to be expressed in mm
    // while the energy in MeV
    primaryGenerator->setParticle(name);
    primaryGenerator->setMomentum(dir);
    primaryGenerator->setPosition(p);
    primaryGenerator->setEnergy(ke);
 
    // Run geant4
    m_runManager->BeamOn();  

    return StatusCode::SUCCESS;
}

//------------------------------------------------------------------------------
StatusCode G4Generator::finalize() 
{
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "finalize: " << endreq;

    // delete the runManager of geant4
    delete m_runManager;
    
    return StatusCode::SUCCESS;
}




