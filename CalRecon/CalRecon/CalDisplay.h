
#ifndef __CALDISPLAY_H
#define __CALDISPLAY_H 1

#include "GaudiKernel/Algorithm.h"
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"

/**
    Display the Cal reconstructed information
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


private:
    
	IGlastDetSvc* detSvc;


};




#endif