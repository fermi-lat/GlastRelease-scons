// $Header$


// includes
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/ToolFactory.h"
#include "gaudiKernel/Property.h"


#include "GuiSvc/IGuiTool.h"

#include "gui/SubMenu.h"
#include "gui/GuiMgr.h"
#include "gui/Command.h"
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"

#include "World.h"
#include "DisplayGeometry.h"
#if 0
#include "MaterialDisplay.h"
#endif

/** 
* @class DetectorDisplay
*
* @brief  A Tool that sets up the detector display 
* $Header$
*/
class DetectorDisplay : public AlgTool, virtual public IGuiTool 
{
public:

    //! Create the tool
    //! @param name The name of the service
    //!
    DetectorDisplay ( const std::string& type, const std::string& name, const IInterface* parent);

    virtual ~DetectorDisplay (){};

    /// implement to define gui elements: will be called from GuiSvc
    StatusCode initialize(gui::GuiMgr*);

private:
    // following specifically for the detector display
    FloatProperty m_axisSize;   //! axis size to use (in mm) 0 for no axis display
    IntegerProperty m_showLevel;    //! number of show commands to reveal details
    DisplayGeometry * m_glast;
};


// declare the service factories for the DetectorDisplay
static ToolFactory<DetectorDisplay> a_factory;
const IToolFactory& DetectorDisplayFactory = a_factory;

/// Standard Constructor
DetectorDisplay::DetectorDisplay(const std::string& type, 
                                 const std::string& name, 
                                 const IInterface* parent)
                                 : AlgTool( type, name, parent ) 
{

    // declare the properties
    declareProperty("axisSize", m_axisSize = 1000);
    declareProperty("showLevel", m_showLevel = 2);

    // Declare additional interface
    declareInterface<IGuiTool>(this);
}


// initialize
StatusCode DetectorDisplay::initialize (gui::GuiMgr* guiMgr) 
{
    StatusCode  status = StatusCode::SUCCESS;

    // bind all of the properties for this service
    status = setProperties ();

    // open the message log
    MsgStream log( msgSvc(), name() );


    // Get the Glast detector service, needed for constants, to keep the GlastDetector pointer
    IGlastDetSvc* gsv;
    if( service( "GlastDetSvc", gsv).isFailure() ) {
        log << MSG::ERROR << "Couldn't set up GlastDetSvc!" << endreq;
        return StatusCode::FAILURE;
    }

    // create the world instance, 
    World* world = World::instance();

    log << MSG::DEBUG << "Instantiating the geometry using the GlastDetSvc " << endreq;
    m_glast = new DisplayGeometry();
    gsv->accept(*m_glast);
    if( log.level() <= MSG::DEBUG ) {
        log << MSG::INFO; 
        m_glast->printStats(log.stream());
        log << endreq;
    }

    // use the display from the Gui service to add the reps
    guiMgr->display().setAxisSize(m_axisSize);



    gui::DisplayRep * worldRep= world->detectorViewer();
    guiMgr->display().add( worldRep, "Detector",  -1 );
    worldRep->update();  // necessary to finish building nested reps
    for(int i =0; i< m_showLevel; ++i) worldRep->show(false); 
#if 0 // did not port: needs some logic to save material names
    // special display of all Mediums with a given material
    MaterialDisplay mdxxx( guiMgr );
#endif
    return status;
}



