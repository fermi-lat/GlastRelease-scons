#ifndef IAcdCalibSvc_H
#define IAcdCalibSvc_H
//  $Header$

// LOCAL INCLUDES

// GLAST INCLUDES
#include "CalibData/Acd/AcdPed.h" 
#include "CalibData/Acd/AcdGain.h"
#include "idents/AcdId.h" 

// EXTLIB INCLUDES
#include "GaudiKernel/IInterface.h"

// STD INCLUDES
//#include <vector>

// Declaration of the interface ID ( interface id, major version,
// minor version)
static const InterfaceID IID_IAcdCalibSvc("IAcdCalibSvc", 1, 0);

/*! @class IAcdCalibSvc
 * \brief Abstract interface for provision of GLAST LAT ACD calibration constants
 * \author Eric Charles (from Zach Fewtrell's CalCalibSvc)
 *
 * \note functions are provided for each calibration type.  
 *       calib objects are passed back by reference to pointer
 *
 */

namespace AcdUtil {
  
  class IAcdCalibSvc : virtual public IInterface {
  public:
    static const InterfaceID& interfaceID() { return IID_IAcdCalibSvc; }
    

    /** \brief get pedestals for given channel
	\param id  the tile or ribbon id
	\param pmt A(0) or B(1) pmt
	\param pedestal a pointer to the relevent pedestal data
    */
    virtual StatusCode getPedestal(idents::AcdId id, unsigned pmt,
				   CalibData::AcdPed*& pedestal) = 0;
			      
    /** \brief get mip peak for a given channel
	\param id  the tile or ribbon id
	\param pmt A(0) or B(1) pmt
	\param pedestal a pointer to the relevent gain data
    */
    virtual StatusCode getMipPeak(idents::AcdId id, unsigned pmt,
				  CalibData::AcdGain*& mipPeak) = 0;

  };
};

#endif // IAcdCalibSvc_H
