
#ifndef __CalValsCorrTool_H
#define __CalValsCorrTool_H 1

#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "EnergyCorr.h"

class IValsTool;

/**   
* @class CalValsCorrTool
*
* Algorithm for reconstruction of energy and direction of incident particle
*
*
* Performs high level energy corrections
*
* The reconstruction here uses CalXtalRecCol to produce a CalClusterCol.
*  It evaluates the barycenter for each layer using the coordinates stored in
*  the CalXtalRecCol,
*
* $Header$
*/


class CalValsCorrTool : public AlgTool, public EnergyCorr {

public:
    
    //! destructor
    CalValsCorrTool( const std::string& type, const std::string& name, const IInterface* parent);
     ~CalValsCorrTool() {}; 
    
     StatusCode initialize();

        
           
     StatusCode doEnergyCorr(double eTotal, Event::CalCluster* cluster);

     StatusCode finalize();
    
     StatusCode execute();

    
private:

    /// name of Tool for finding last layer energy leakage
    std::string m_calValsToolName;

    /// pointer to actual tool for last layer energy correlation
    IValsTool* m_calValsTool;



};

#endif



