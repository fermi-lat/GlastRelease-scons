
#ifndef __CALDISPLAY_H
#define __CALDISPLAY_H 1

#include "GaudiKernel/Algorithm.h"

class CalRecLogs;
class CsIClusterList;
class ICalGeometrySvc;

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

    void clearLogsDisp(){m_crl=0;}
    void clearClusterDisp(){m_cls=0;}

private:
    
	CalRecLogs* m_crl;
	CsIClusterList* m_cls;
	ICalGeometrySvc* m_CalGeo; 

};




#endif