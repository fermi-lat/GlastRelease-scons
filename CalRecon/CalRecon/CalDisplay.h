
#ifndef __CALDISPLAY_H
#define __CALDISPLAY_H 1

#include "GaudiKernel/Algorithm.h"
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "gui/DisplayControl.h"

/**
* @class CalDisplay
*  
* @brief This class provides the display of the Cal reconstructed information
*
* CalDisplay algorithm displays the results of calorimeter reconstruction
* from TDS classes CalXtalRecCol and CalClustersCol. See details in
*  CalRep class description.   
* Display is set up by initialize() method, which adds CalRep class object
* into the list of objects to be displayed.
*
*
* @author A.Chekhtman
*
* $Header$
*/
class CalDisplay : public Algorithm
{
public:

    /// constructor required for Gaudi algorithm
    CalDisplay(const std::string& name, ISvcLocator* pSvcLocator);
    
    virtual ~CalDisplay() {}
    
    /**
    *   This function does the initialization of the calorimeter display:
    *
    *    - extracts necessary detector geometry constants using GlastDetSvc;
    * 
    *    - gets information on Z position of top and bottom calorimeter
    *      layers;
    *
    *    - creates the object of CalRep class and, using GiuSvc, adds
    *      the pointer to this object to the listto be displayed by
    *      graphic user interface. 
    *       
    */

    StatusCode initialize();

    /// empty function, required for Gaudi algorithm    
    StatusCode execute(){ return StatusCode::SUCCESS;}

    /// empty function, required for Gaudi algorithm    
    StatusCode finalize(){ return StatusCode::SUCCESS;}


private:
    
    /// pointer to GlastDetSvc, used to get detector geometry constants
    IGlastDetSvc* detSvc;
};

/** 
*  @class CalRep
*
*  @brief This class provides graphic representation of calorimeter
*         reconstructed data.
*
*  It draws following elements on gui display
*
*  - red box in each crystal with non-zero energy deposition,
*    box size is propotional to the energy in this crystal, box center
*    shows the position reconstructed from signal asymmetry;
*
*  - blue cross in the average position calculated for each layer;
*
*  - green diamond in the average position calculated for the cluster;
*
*  - green line showing the reconstructed shower direction;  
*
*  @author A.Chekhtman
*/
class CalRep : public gui::DisplayRep {
    
public:

    /// constructor with parameters to initialize private data members    
    CalRep(IDataProviderSvc* eventSvc, float xtalHeight,
           float calZtop, float calZbottom)
        :m_eventSvc(eventSvc),m_xtalHeight(xtalHeight),
         m_calZtop(calZtop), m_calZbottom(calZbottom){}


    /// function making the real drawing,
    /// it is called when display is updated
    void update();

private:

    /// pointer to EventSvc
    IDataProviderSvc* m_eventSvc;

    /// crystal height
    float m_xtalHeight;

    /// Z position of top calorimeter layer
    float m_calZtop;

    /// Z position of bottom calorimeter layer
    float m_calZbottom;
};
#endif