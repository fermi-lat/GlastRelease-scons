// File and Version Information:
// $Header$
// Description:
// Handles the display of the ACD DOCA values on the GUI.

#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"

#include "Event/TopLevel/EventModel.h"
#include "Event/Recon/AcdRecon/AcdRecon.h"
#include "idents/VolumeIdentifier.h"

#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "CLHEP/Geometry/Point3D.h"

// gui, display includes
#include "GuiSvc/IGuiSvc.h"
#include "gui/DisplayControl.h"
#include "gui/GuiMgr.h"

/** @class AcdDisplay
 * @brief Display the ACD DOCA's as numbers at the tile center positions
 *
 * @author Heather Kelly
 * $Header$
 */
class AcdDisplay : public Algorithm
{
public:
    AcdDisplay(const std::string& name, ISvcLocator* pSvcLocator); 
    virtual ~AcdDisplay() {}

    StatusCode initialize();
    StatusCode execute(){ return StatusCode::SUCCESS;}
    StatusCode finalize(){ return StatusCode::SUCCESS;}
    class DocaRep;
};

static const AlgFactory<AcdDisplay>  Factory;
const IAlgFactory& AcdDisplayFactory = Factory;

AcdDisplay::AcdDisplay(const std::string& name, ISvcLocator* pSvcLocator) :
Algorithm(name, pSvcLocator) 
{
    
}

/** @class AcdRep
 * @brief Handles the mechanics of updating the display
 *
 * @author Heather Kelly
 */
class AcdDisplay::DocaRep : public gui::DisplayRep {
public:
    DocaRep(IDataProviderSvc* dps, IGlastDetSvc* detsvc):m_dps(dps), m_detsvc(detsvc){}
    void update(){
        SmartDataPtr<Event::AcdRecon> acdRec(m_dps, EventModel::AcdRecon::Event);
        if (!acdRec) return;

	// loop thru the tiles with DOCA's and show marker and value
        // no, only show the closest one for now
	if( acdRec->getDoca()>1000. ) return;

        // get the VolumeIdentifier of the closest doca (TODO: move this to idents::AcdId)
        idents::AcdId id = acdRec->getMinDocaId();
        idents::VolumeIdentifier volid;
        volid.append(1); 
        volid.append(id.face()); 
        volid.append(id.column()); 
        volid.append(id.row());

        HepPoint3D tilecenter;
        HepTransform3D transform;

        StatusCode sc = m_detsvc->getTransform3DByID(volid, &transform);
        if(sc.isSuccess()) tilecenter = transform*tilecenter;

	set_color("red");
	markerAt(tilecenter);
	move_to(tilecenter); 
	char text[10]; sprintf(text, "%8g", acdRec->getDoca());
	drawText(text);
	set_color("black");
    }
private:
    IDataProviderSvc* m_dps;
    IGlastDetSvc* m_detsvc;
};

StatusCode AcdDisplay::initialize()
{
    //Look for the gui service
    IGuiSvc*   guiSvc = 0;
    StatusCode sc     = service("GuiSvc", guiSvc);
    if( sc.isFailure() )  return sc;

    // need the detector service to find position of tiles by id.
    IGlastDetSvc* detsvc = 0;
    if( service( "GlastDetSvc", detsvc).isFailure() ) return StatusCode::FAILURE;
  
    gui::DisplayControl& display = guiSvc->guiMgr()->display();
    
    gui::DisplayControl::DisplaySubMenu& acdmenu = display.subMenu("AcdRecon");
    
    acdmenu.add(new DocaRep(eventSvc(), detsvc), "doca");
    // TODO: add other reps as needed.
    
    return StatusCode::SUCCESS;
}

