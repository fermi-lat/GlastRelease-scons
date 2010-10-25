#ifndef ICalClassifyTool_h
#define ICalClassifyTool_h 

#include "GaudiKernel/IAlgTool.h"
#include "Event/Recon/CalRecon/CalCluster.h"
#include <vector>

/**   
* @class ICalClassifyTool
*
* Base class for CAL cluster classifier tools
*
* $Header: 
*/

static const InterfaceID IID_ICalClassifyTool("ICalClassifyTool",1,0) ;
    
//! Typedefs to define associated Xtal's in cluster analysis -- TBD something ?
//typedef  std::list<Event::CalXtalRecData *> XtalDataList ;
//typedef  std::vector<XtalDataList *>        XtalDataListVec ;

class ICalClassifyTool : virtual public IAlgTool {

  public:

    // retrieve Gaudi interface ID
    static const InterfaceID& interfaceID() {return IID_ICalClassifyTool;}

    //! main method
    virtual StatusCode classifyClusters(Event::CalClusterCol* calClusterCol) = 0;
 } ;

#endif



