
#ifndef __CalSingleClusteringTool_H
#define __CalSingleClusteringTool_H 1

#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "CalClusteringTool.h"

/**   
* @class CalSingleClusteringTool
*
* Find single cluster from all CAL hits
*
*
* $Header$
*/


class CalSingleClusteringTool : public CalClusteringTool
 {
  public:
    
    CalSingleClusteringTool
     ( const std::string & type, 
       const std::string & name, 
       const IInterface * parent ) ;
    ~CalSingleClusteringTool() ;
        
  protected:

    //! Distinguish sets of related xtals
    virtual void makeSets( const XtalDataVec & xtals, XtalDataVecVec & clusters ) ;

 } ;

#endif



