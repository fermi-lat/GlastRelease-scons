#ifndef IXtalDigiTool_H
#define IXtalDigiTool_H
// $Header$
/*! @class IXtalDigiTool
 *
 * \author Zach Fewtrell
 *
 * \brief Abstract interface class for estimation of the digital response of 
 one calorimeter crystal to a set of energy depositions.
 * 
 *
 */

// GLAST INCLUDES
#include "Event/MonteCarlo/McIntegratingHit.h"
#include "Event/TopLevel/Event.h"
#include "Event/Digi/CalDigi.h"
#include "Event/Digi/GltDigi.h"
#include "CalUtil/CalDefs.h"
#include "CalUtil/CalArray.h"

// EXTLIB INCLUDES
#include "GaudiKernel/IAlgTool.h"

// STD INCLUDES
#include <vector>

static const InterfaceID IID_IXtalDigiTool("IXtalDigiTool", 1, 0);

using namespace std;
using namespace idents;
using namespace CalUtil;


class IXtalDigiTool : virtual public IAlgTool {
 public:

  static const InterfaceID& interfaceID() { return IID_IXtalDigiTool; }

  /** \brief calculate Adc response for one cal xtalIdx.  also select best rng.

  \param hitList input vector of energy depositions.  
  \param evtHdr optional pointer to current event header (used for runid & evtid)
  \param calDigi output empty CalDigi object to be populated (xtalId & rangeMode expected populated)
  \param lacBits output boolean for log accept on each xtal face
  \param trigBits output booleans for triggers for each xtal diode.
  \param glt optional output GltDigi class.  Will populate if (glt != 0).

  */
  virtual StatusCode calculate(const vector<const Event::McIntegratingHit*> &hitList,
                               const Event::EventHeader *evtHdr,            
                               Event::CalDigi &calDigi,
                               CalArray<FaceNum, bool> &lacBits,
                               CalArray<XtalDiode, bool> &trigBits,
                               Event::GltDigi *glt
                               ) = 0;

};

#endif //_IXtalDigiTool_H
