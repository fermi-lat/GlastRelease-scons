#ifndef ICalTrigTool_H
#define ICalTrigTool_H
// $Header$
/*! @class ICalTrigTool
 *
 * \author Zach Fewtrell
 *
 * \brief Abstract interface class for estimation of the FLE & FHE trigger 
 response of GLAST calorimeter w/ given digi ADC response.
 *
 * provides functions for both single xtal & entire Cal
 * 
 * \note Incomplete/single range digi readouts will prevent simulation of 
 certain rare - direct diode deposit related results such as a high FHE trigger 
 on a crystal w/ low range ULD readout.
 * 
 *
 */

// GLAST INCLUDES
#include "Event/TopLevel/Event.h"
#include "CalUtil/CalDefs.h"
#include "CalUtil/CalArray.h"
#include "Event/Digi/CalDigi.h"
#include "Event/Digi/GltDigi.h"

// EXTLIB INCLUDES
#include "GaudiKernel/IAlgTool.h"
#include "GaudiKernel/IDataProviderSvc.h"

// STD INCLUDES
#include <vector>

static const InterfaceID IID_ICalTrigTool("ICalTrigTool", 1, 0);

using namespace std;
using namespace idents;
using namespace CalUtil;

class ICalTrigTool : virtual public IAlgTool {
 public:

  static const InterfaceID& interfaceID() { return IID_ICalTrigTool; }


  /** \brief calc Trigger response for single cal xtal, given 4 range adc output

  \param xtalIdx specify Cal xtal log
  \param adcPed input xtal digi readout 1 ped subtracted adc per xtal face/rng 
  \param trigBits output trigger bits.. 1 per xtal diode.
  \param glt optional output GltDigi class.  Will populate if (glt != 0).

  */
  virtual StatusCode calcXtalTrig(XtalIdx xtalIdx,
                                  const CalArray<XtalRng, float> &adcPed,
                                  CalArray<XtalDiode, bool> &trigBits,
                                  Event::GltDigi *glt
                                  ) = 0;

  /** \brief calc Trigger response for single cal xtal
  
  \param calDigi input ADC readout info
  \param trigBits output trigger bits.. 1 per xtal diode.
  \param glt optional output GltDigi class.  Will populate if (glt != 0).

  \note Incomplete/single range digi readouts will prevent simulation of 
  certain rare - direct diode deposit related results such as a high FHE 
  trigger on a crystal w/ low range ULD readout.
  */
  virtual StatusCode calcXtalTrig(const Event::CalDigi& calDigi,
                                  CalArray<XtalDiode, bool> &trigBits,
                                  Event::GltDigi *glt
                                  ) = 0;


  /** \brief calc Trigger response for sigle Cal xtal readout
      
  \param xtalIdx specify Cal xtal log
  \param ro single Cal Xtal ADC readout
  \param trigBits output trigger bits.. 1 per xtal diode.
  \param glt optional output GltDigi class.  Will populate if (glt != 0).

  \note Incomplete/single range digi readouts will prevent simulation of 
  certain rare - direct diode deposit related results such as a high FHE 
  trigger on a crystal w/ low range ULD readout.
  */
  virtual StatusCode calcXtalTrig(XtalIdx xtalIdx,
                                  const Event::CalDigi::CalXtalReadout &ro,
                                  CalArray<XtalDiode, bool> &trigBits,
                                  Event::GltDigi *glt
                                  ) = 0;

  /** \brief global calc CAL FLE & FHE trigger response based current
      Event Digi info

      \param calDigiCol input collection of CalDigi objects
      \param glt optional output GltDigi class.  Will populate if (glt != 0).
      \param fle output Cal global FLE trigger bit
      \param fhe output Cal global FHE trigger bit

      \note Incomplete/single range digi readouts will prevent simulation of 
      certain rare - direct diode deposit related results such as a high FHE 
      trigger on a crystal w/ low range ULD readout.
  */
  virtual StatusCode calcGlobalTrig(const Event::CalDigiCol& calDigiCol,
                                    CalArray<DiodeNum, bool> &trigBits,
                                    Event::GltDigi *glt
                                    ) = 0;

  /** \brief search for GltDigi class in TDS & create a new one if needed.
   */
  virtual Event::GltDigi* setupGltDigi(IDataProviderSvc *eventSvc) = 0;


};

#endif //
