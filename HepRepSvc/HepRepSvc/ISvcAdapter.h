// $Header:

#ifndef ISVCADAPTER_H
#define ISVCADAPTER_H

#include <vector>
#include <string>

/* @class: ISvcAdapter 
 *
 * @brief: The abstract interface to the Gaudi service adapter
 *
 * @author:R.Giannitrapni
 */ 

class ISvcAdapter
{
 public:
  /// used to step to the next event 
  virtual void nextEvent(int) = 0;
  /// return a list of sources names to be used by FluxSvc
  virtual std::vector<std::string>& getSources() = 0;
  /// set the source to be used by FluxSvc
  virtual void setSource(std::string) = 0;
};

#endif //ISVCADAPTER_H
