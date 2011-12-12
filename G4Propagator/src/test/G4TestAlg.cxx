// $Header$

// Include files

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/GaudiException.h" 
#include <string>

// This defines the interface for the propogator
#include "GlastSvc/Reco/IPropagator.h"

// This gets the GLAST detector service for access to geometry
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "idents/VolumeIdentifier.h"
#include "idents/TowerId.h"
#include "CLHEP/Geometry/Transform3D.h"


/*! \class G4TestAlg
\brief 
  This is a place to put test code to examine what G4Generator has done.
  */

class G4TestAlg : public Algorithm {
    
public:
    //! Constructor of this form must be provided
    G4TestAlg(const std::string& name, ISvcLocator* pSvcLocator); 
    
    StatusCode initialize();
    StatusCode execute();
    StatusCode finalize();
    
private:
    IPropagator*   m_propagator;    

    /// pointer to the detector service
    IGlastDetSvc * m_pDetSvc;

    double         m_startx;
    double         m_starty;
};


//static const AlgFactory<G4TestAlg>  Factory;
//const IAlgFactory& G4TestAlgFactory = Factory;
DECLARE_ALGORITHM_FACTORY(G4TestAlg);

//
G4TestAlg::G4TestAlg(const std::string& name, ISvcLocator* pSvcLocator) :
Algorithm(name, pSvcLocator)
{
    declareProperty("StartX", m_startx = 200.); 
    declareProperty("StartX", m_starty = 200.); 
}


/*! */
StatusCode G4TestAlg::initialize() {
      

    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "initializing..." << endreq;
    StatusCode  sc = StatusCode::SUCCESS;
    
    // Use the Job options service to set the Algorithm's parameters
    setProperties();

    //Locate a pointer to the G4Propagator
    //IPropagator* propagatorTool = 0;
    if( (sc = toolSvc()->retrieveTool("G4PropagationTool", m_propagator)).isFailure() )
    {
        throw GaudiException("ToolSvc could not find G4PropagationTool", name(), sc);
    }

    // Locate the GlastDetSvc
    if ( (sc = service("GlastDetSvc", m_pDetSvc)).isFailure() )
    {
        throw GaudiException("Could not find GlastDetSvc", name(), sc);
    }

    double stayClear = 0.;
    if( (sc = m_pDetSvc->getNumericConstByName("TKRVertStayClear", &stayClear)).isFailure()) 
    {
        throw GaudiException("Couldn't find TKRVertStayClear", name(), sc);
    }

    // get the z coordinate of the top silicon in the bottom tray 
    //         (should be pretty safe!)
    // don't make any assumptions about the view in the bottom tray
    idents::VolumeIdentifier bottom;
    
    bottom.init(0,0);
    bottom.append(0);        // in Tower
    idents::TowerId t(1);    //
    bottom.append(t.iy());   // yTower
    bottom.append(t.ix());   // xTower
    bottom.append(1);        // Tracker
    bottom.append(0);        // tray 0
    idents::VolumeIdentifier idBot;
    HepTransform3D botTransform;
    for (int view = 0; view<2; ++view) 
    {
        idBot = bottom;
        idBot.append(view);          // try both views
        idBot.append(1);             // top silicon (*most* bottom trays have one!)
        idBot.append(0); idBot.append(0);  // ladder, wafer
        //std::cout << "view " << view << " idBot " << idBot.name() << std::endl;
        if(sc = m_pDetSvc->getTransform3DByID(idBot, &botTransform).isSuccess()) {
            break;
        }
    }

    // Everything ok?
    if (sc.isFailure()) 
    {
        throw("Couldn't find bottom tray silicon", name(), sc);
    }

    double zBot;
    zBot = (botTransform.getTranslation()).z();

    // this is always above the tracker
    double propTop = stayClear+zBot;

    // Define starting point and direction
    Point  startPoint(m_startx, m_starty, propTop);
    Vector startDir(0., 0., -1.);

    // Embed call to propagator in a try-catch block 
    try
    {
        // Initialize the propagator with starting point and direction
        // This causes the propagator to locate itself in the G4 geometry
        m_propagator->setStepStart(startPoint, startDir);

        // Determine an arclength that will run the length of the tracker in z
        double propRange = stayClear+100.;
        double propBot = propTop - propRange;

        // Call the propagator to now step from the starting position through the given
        // arclength in the defined direction
        m_propagator->step(propRange);

        log << MSG::INFO  << "Propagator goes from "<< propTop << " to " << propBot << endreq;

        int numSteps = m_propagator->getNumberSteps();

        log << MSG::INFO  << "Propagator took " << numSteps << " steps" << endreq;
        //int istep;
        //idents::VolumeIdentifier id;
        //idents::VolumeIdentifier prefix = m_pDetSvc->getIDPrefix();

        // This will print out the results of each step along the way
        log << MSG::INFO;
        m_propagator->printOn(log.stream());
        log << endreq;
    }
    catch(...)
    {
        log << MSG::INFO << "G4TestAlg has caught an exception thrown during propagation, passing along... " << endreq;
        throw;
    }

    return sc;
}


StatusCode G4TestAlg::execute() {
    
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream   log( msgSvc(), name() );

    return sc;
}


StatusCode G4TestAlg::finalize() {
    
    return StatusCode::SUCCESS;
}






