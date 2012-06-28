#ifndef CalClusterMap_H
#define CalClusterMap_H


#include <string>
#include <map>
#include "Event/Recon/CalRecon/CalClusterVec.h"
#include "Event/TopLevel/EventModel.h"


static const CLID& CLID_CalClusterMap = InterfaceID("CalClusterMap",  1, 0);

/** 
* @class CalClusterMap
*
* @brief Gaudi TDS class providing access to the CAL clusters for the event
* level analysis.
*
* This is essentially a std::map of std::vectors of pointers to
* Event::CalClusterVec objects. While we still keep an Event::CalClusterCol a
* the basic container for the cluster objects (i.e. owning the objects
* themselves and taking care of the related cleanup) we won't use that
* container for accessing the actual cluster, this new class providing
* a more flexible access option.
*
* We decided to map vectors of clusters rather than clusters as a more general
* data structure, allowing to group together in a straighforward way, e.g., the
* "raw" output of the clustering stage.
*
* @author Luca Baldini.
*/


namespace Event { //Namespace Event
 
  class CalClusterMap : virtual public std::map<std::string, CalClusterVec>,
    virtual public DataObject
    {
    public:
      /// Constructor.
    CalClusterMap() : DataObject() {clear();}

      /// Destructor.
      virtual ~CalClusterMap() {};

      /// Some more Gaudi-related stuff.
      virtual const CLID& clID() const  { return CalClusterMap::classID(); }
      static const CLID& classID()      { return CLID_CalClusterMap; }

      /// Return the CalClusterVec object corresponding to a given key.
      CalClusterVec get(std::string);
      /// Return the CalClusterVec corresponding to the raw output from the clustering.
      CalClusterVec getRawClusterVec() { return get(EventModel::CalRecon::CalRawClusterVec); }
      /// Return a pointer to the cluster at the front of a given CalClusterVec.
      /// This is a convenience function for accessing the vectors with a single element
      /// (e.g., the one containing the user cluster).
      CalCluster* getFront(std::string key) { return get(key).front(); }      
      /// Return a pointer to the uber cluster.
      CalCluster* getUberCluster() { return getFront(EventModel::CalRecon::CalUberCluster); }
  };
  
  // Iterator typedefs.
  typedef CalClusterMap::iterator       CalClusterMapItr;
  typedef CalClusterMap::const_iterator CalClusterMapConItr;


}; //Namespace Event

#endif        
