/**
 * @class ICalClusterHitTool
 *
 * @brief Implements an interface for a Gaudi Tool for having access to the
 *        calorimeter hit xtals associated with a given cluster (based on the)
 *        information in the TDS relational tables.
 *
 * @author Luca Baldini (luca.baldini@pi.infn.it)
 *
 * $Header$
 */

#ifndef ICalClusterHitTool_h
#define ICalClusterHitTool_h

#include "GaudiKernel/IAlgTool.h"
#include "Event/Recon/CalRecon/CalCluster.h"

static const InterfaceID IID_ICalClusterHitTool("ICalClusterHitTool", 1 , 0);

class ICalClusterHitTool : virtual public IAlgTool
{

 public:
  
  /// Define the interfaces for the derived classes.
  virtual StatusCode fillRecDataVec(Event::CalCluster*) = 0;
  virtual std::vector<Event::CalXtalRecData*> getRecDataVec() const = 0;
  virtual void setRecDataVec(std::vector<Event::CalXtalRecData*>) = 0;

  /// Retrieve interface ID
  static const InterfaceID& interfaceID() { return IID_ICalClusterHitTool; }
};

#endif
