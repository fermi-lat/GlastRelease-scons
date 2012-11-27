#ifndef IXtalRecTool_H
#define IXtalRecTool_H
// $Header$
/** @file
    @author Z.Fewtrell
*/
/*! @class IXtalRecTool
 *
 * \author Z.Fewtrell
 *
 * \brief Abstract interface class for calculation of total deposited energy 
 & centroid position for a single calorimeter crystal.
 * 
 *
 */

// LOCAL INCLUDES

// GLAST INCLUDES
#include "CalUtil/CalVec.h"
#include "CalUtil/CalDefs.h"

// EXTLIB INCLUDES
#include "GaudiKernel/IAlgTool.h"

// STD INCLUDES

static const InterfaceID IID_IXtalRecTool("IXtalRecTool", 1, 3);

// forward declarations
namespace Event {
	class CalDigi;
	class CalXtalRecData;
};

class INeighborXtalkTool;

class IXtalRecTool : virtual public IAlgTool {

 public:
  static const InterfaceID& interfaceID() { return IID_IXtalRecTool; }

  /** \brief estimate total deposited energy given the digital
      response for both faces

      \param digi single input digi readout 

      \param xtalRec expected to contain calxtalid field but no 
      CalRangeRecData objects.  calculate() will attempt 
      to create 1 rangeRec object, but this will not always be possible.
      Caller should check xtalRec.getNReadouts() to see if xtal recon 
      was successful.  It _is_ possible for this function to return 
      StatusCode::SUCCESS along w/ 0 CalRangeRecData objects.

      \param belowThresh returns TRUE for each face signal below noise for any 
      adc range

      \param xtalBelowThresh returns TRUE if either face signal below LEX8 noise level

      \param saturated returns TRUE for each face signal that is >= HEX1 
      saturation level

	  \param nbrXtalkTool (optional) pointer to neighboring xtal->xtal electronic crosstalk model

      \return It _is_ possible for this function to return 
      StatusCode::SUCCESS along w/ 0 CalRangeRecData objects.


  */
  virtual StatusCode calculate(const Event::CalDigi &digi,
                               Event::CalXtalRecData &xtalRec,
                               CalUtil::CalVec<CalUtil::FaceNum, bool> &belowNoise,
                               CalUtil::CalVec<CalUtil::FaceNum, bool> &saturated,
                               CalUtil::CalVec<CalUtil::FaceNum, bool> &adcSaturated,
                               INeighborXtalkTool const*const nbrXtalkTool
                               ) = 0;
};

#endif //_IXtalRecTool_H
