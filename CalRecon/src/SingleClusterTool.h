
#ifndef __SingleClusterTool_H
#define __SingleClusterTool_H 1

#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "Cluster.h"

/**   
* @class SingleClusterTool
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



