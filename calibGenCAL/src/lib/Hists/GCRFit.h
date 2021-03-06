#ifndef GCRFit_h
#define GCRFit_h

// $Header$

/** @file
    @author Zachary Fewtrell
*/

// LOCAL INCLUDES

// GLAST INCLUDES

// EXTLIB INCLUDES

// STD INCLUDES
#include <string>

class TDirectory;

namespace CalUtil {
  class CalMPD;
}

namespace calibGenCAL {
  class GCRHists;

  /** \brief Collection of tools for fitting GCR calib histograms
   */
  namespace GCRFit {
    /// fit GCR histograms & output fitting results to CalMPD obj w/ simple Gaussian peakshape
    /// \parm calMPD output calibration constants
    /// \parm writeFile location for output fit results tuple
    /// \parm tupleName output tuple fit results name
    void gcrFitGaus(GCRHists &histCol,
			CalUtil::CalMPD &calMPD,
                    TDirectory *const writeFile,
                    const std::string &tupleName="GCRFitGauss"
                    );
  } // namespace GCRFit

} // namespace calibGenCAL 

#endif
