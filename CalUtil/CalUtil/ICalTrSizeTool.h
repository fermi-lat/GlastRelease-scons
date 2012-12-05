/**
 * @class ICalDirFittingTool
 *
 * @brief Implements a Gaudi Tool for evaluating the transverse size of
 *        a CAL cluster, give an input axis.
 *
 * @author Luca Baldini (luca.baldini@pi.infn.it)
 *
 * $Header$
 */

#ifndef ICalTrSizeTool_h
#define ICalTrSizeTool_h

#include "GaudiKernel/IAlgTool.h"
#include "Event/Recon/CalRecon/CalXtalRecData.h"
#include "Event/Recon/CalRecon/CalCluster.h"


static const InterfaceID IID_ICalTrSizeTool("ICalTrSizeTool", 1 , 0);

class ICalTrSizeTool : virtual public IAlgTool
{
 public:
  /// Define the interfaces for the derived classes.
  virtual StatusCode fill(std::vector<Event::CalXtalRecData*> xtalList) = 0;
  virtual StatusCode fill(Event::CalCluster* cluster) = 0;
  virtual StatusCode computeTrSize(const Point& origin,
                                   const Vector& direction) = 0;
  virtual StatusCode computeTrSizeTrans(const Point& origin,
                                        const Vector& direction) = 0;
  virtual double getQuantile(double frac) = 0;

  /// Retrieve interface ID
  static const InterfaceID& interfaceID() { return IID_ICalTrSizeTool; }

 private:
  virtual StatusCode compute(const Point& origin, const Vector& direction,
                             bool transverse) = 0;
};

#endif
