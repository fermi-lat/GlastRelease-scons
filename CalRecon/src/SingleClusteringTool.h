
#ifndef __SingleClusteringTool_H
#define __SingleClusteringTool_H 1

#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "ClusteringTool.h"

/**   
* @class SingleClusteringTool
*
* Find single cluster from all CAL hits
*
*
* $Header$
*/


class SingleClusteringTool : public ClusteringTool
 {
  public:
    
    SingleClusteringTool
     ( const std::string & type, 
       const std::string & name, 
       const IInterface * parent ) ;
    ~SingleClusteringTool() ;
        
  protected:

    // use a single global cluster
    xTalDataVec nextXtalsSet( xTalDataVec & xTalVec ) ;

 } ;

#endif



