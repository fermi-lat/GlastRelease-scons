
#ifndef __SingleClusterTool_H
#define __SingleClusterTool_H 1

#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "Cluster.h"

/**   
* @class SingleClusterTool
*
* Find single cluster from all CAL hits
*
*
* $Header$
*/


class SingleClusterTool : public AlgTool, public Cluster {

public:
    
    //! destructor
    SingleClusterTool( const std::string& type, const std::string& name, const IInterface* parent);
     ~SingleClusterTool() {}; 
    
     StatusCode initialize();

        
/*!Performs the reconstruction, creates one CalCluster object and stores
 * there the following results: 
 * - Energy per layer is computed and stored in CalCluster in MeV
 * - Barycenter per layer is also computed and stored in CalCluster
 */        
     StatusCode findClusters(Event::CalXtalRecCol* calXtalRecCol);

     StatusCode finalize();
    
     StatusCode execute();

protected:
     // calculate the shower direction from the CAL hits
    Vector Fit_Direction(std::vector<Vector> pos,
                                     std::vector<Vector> sigma2,
                                     int nlayers);
    
    
    
private:
       
	//! crystal width
    double m_CsIWidth;

    //! crystal height
    double m_CsIHeight;

    //! number of layers
    int m_CalnLayers;

    //! pointer to GlasDetSvc
    IGlastDetSvc* detSvc;

};

#endif



