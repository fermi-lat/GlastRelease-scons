#ifndef _IXtalRecTool_H
#define _IXtalRecTool_H 1
/*! @class IXtalRecTool
 *
 * \author Zach Fewtrell
 *
 * \brief Abstract interface class for calculation of total deposited energy 
 & centroid position for a single calorimeter crystal.
 * 
 *
 */

// LOCAL INCLUDES
#include "CalXtalResponse/CalTupleEntry.h"


// GLAST INCLUDES
#include "idents/CalXtalId.h"
#include "Event/TopLevel/Event.h"
#include "Event/Digi/CalDigi.h"
#include "Event/Recon/CalRecon/CalXtalRecData.h"

// EXTLIB INCLUDES
#include "GaudiKernel/IAlgTool.h"

// STD INCLUDES

using namespace idents;

static const InterfaceID IID_XtalRecTool("IXtalRecTool", 1 , 0);

class IXtalRecTool : virtual public IAlgTool {
 public:
  static const InterfaceID& interfaceID() { return IID_XtalRecTool; }

  /** \brief estimate total deposited energy given the digital response for both faces
      \param xtalId specify xtal log
      \param readout single input digi readout 
      \param evt pointer to current event (used for runid & evtid)
      \param rngRec output recon readout.
      \param belowThreshP returns TRUE if POS face signal is below noise (for all ranges)
      \param belowThreshN returns TRUE if NEG face signal is below noise (for all ranges)
      \param saturatedP POS face ADC is saturated (for all ranges)
      \param saturatedN NEG face ADC is saturated (for all ranges)
      \param calTupleEnt pointer to optional tuple entry (NULL pointer is ignored)
  */
  virtual StatusCode calculate(const Event::EventHeader &evtHdr,
                               const Event::CalDigi &digi,
                               Event::CalXtalRecData &xtalRec, // output
                               bool &belowThreshP,    // output
                               bool &belowThreshN,    // output
                               bool &xtalBelowThresh,
                               bool &saturatedP,      // output
                               bool &saturatedN,       // output
                               CalTupleEntry *calTupleEnt
                               ) = 0;
};

#endif //_IXtalRecTool_H
