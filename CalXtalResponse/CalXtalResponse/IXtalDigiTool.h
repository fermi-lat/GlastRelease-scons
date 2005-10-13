#ifndef _IXtalDigiTool_H
#define _IXtalDigiTool_H
/*! @class IXtalDigiTool
 *
 * \author Zach Fewtrell
 *
 * \brief Abstract interface class for estimation of the digital response of one calorimeter crystal to a set of energy depositions.
 * 
 *
 */

// GLAST INCLUDES
#include "idents/CalXtalId.h"
#include "Event/MonteCarlo/McIntegratingHit.h"
#include "Event/TopLevel/Event.h"
#include "Event/Digi/CalDigi.h"

// EXTLIB INCLUDES
#include "GaudiKernel/IAlgTool.h"

// STD INCLUDES
#include <vector>

static const InterfaceID IID_XtalDigiTool("IXtalDigiTool", 1 , 0);

using namespace std;
using namespace idents;

class IXtalDigiTool : virtual public IAlgTool {
 public:

  static const InterfaceID& interfaceID() { return IID_XtalDigiTool; }

  /** \brief calculate Adc response for one cal xtalId.  also select best rng.
      \param CalXtalId specify xtal log
      \param hitList input vector of energy depositions.  (const *) is used to save space.
      \param evt pointer to current event (used for runid & evtid)
      \param calDigi output empty CalDigi object to be populated
      \param fleP output FLE trigger Pos face
      \param fleN output FLE trigger Pos face
      \param fheP output FHE trigger Neg face
      \param fheN output FHE trigger Neg face
      \param lacP output boolean for log accept on Positive xtal face
      \param lacN output boolean for log accept on Negative xtal face
  */
  virtual StatusCode calculate(CalXtalId xtalId, 
                               const vector<const Event::McIntegratingHit*> &hitList,
                               const Event::EventHeader &evtHdr,            
                               Event::CalDigi &calDigi,     // output 
                               bool &lacP,                  // output 
                               bool &lacN,                  // output 
                               bool &fleP,                  // output 
                               bool &fleN,                  // output 
                               bool &fheP,                  // output 
                               bool &fheN                   // output 
                               ) = 0;

};

#endif //_IXtalDigiTool_H
