//$Header$

/// Gaudi specific include files
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"

// gui, display includes
#include "GuiSvc/IGuiSvc.h"
#include "gui/DisplayControl.h"
#include "gui/GuiMgr.h"

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/**
    Display the Cal reconstructed centriod, axes
  */
class CalDisplay : public Algorithm
{
public:
    //! Constructor of this form must be provided
    CalDisplay(const std::string& name, ISvcLocator* pSvcLocator); 
    virtual ~CalDisplay() {}
    //! mandatory
    StatusCode initialize();
    //! mandatory
    StatusCode execute();
    //! mandatory
    StatusCode finalize(){ return StatusCode::SUCCESS;}
    
};
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
static const AlgFactory<CalDisplay>  Factory;
const IAlgFactory& CalDisplayFactory = Factory;

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CalDisplay::CalDisplay(const std::string& name, ISvcLocator* pSvcLocator) :
Algorithm(name, pSvcLocator) 
{
    
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class CalRep : public gui::DisplayRep {
public:
    CalRep(){}
    void update(){

		//code goes here
	}
};

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
StatusCode CalDisplay::initialize()
{
    //Look for the gui service
    IGuiSvc* guiSvc = 0;
    StatusCode sc = service("GuiSvc", guiSvc);

    
    //Ok, see if we can set up the display
    if (sc.isSuccess())  {
	//Set up the display rep for Clusters
	guiSvc->guiMgr()->display().add(new CalRep(), "Cal recon");
    }
    
    return sc;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
StatusCode CalDisplay::execute()
{
    StatusCode sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());

    // set pointer to cal data on TDS
    //    ...
    
    return sc;
}


