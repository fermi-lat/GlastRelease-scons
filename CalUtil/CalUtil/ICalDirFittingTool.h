/**
 * @class ICalDirFittingTool
 *
 * @brief Implements a Gaudi Tool for fitting the CAL direction in
 *        several different flavours.
 *
 * @author Luca Baldini (luca.baldini@pi.infn.it)
 *
 * $Header$
 */

#ifndef ICalDirFittingTool_h
#define ICalDirFittingTool_h

#include "GaudiKernel/IAlgTool.h"
#include "Event/Recon/CalRecon/CalCluster.h"

static const InterfaceID IID_ICalDirFittingTool("ICalDirFittingTool", 1 , 0);

class ICalDirFittingTool : virtual public IAlgTool
{

 public:
  
  /// Define the interfaces for the derived classes.
  virtual StatusCode transverseFit2d(Event::CalCluster* cluster,
                                     double powerWeight = 1.) = 0;
  virtual Event::CalFitParams getFitParams() const = 0;

  /// Retrieve interface ID
  static const InterfaceID& interfaceID() { return IID_ICalDirFittingTool; }
};

#endif
