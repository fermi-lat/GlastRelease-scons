
#ifndef __FuzzyClusterTool_H
#define __FuzzyClusterTool_H 1

// The next is the base class of FuzzyCluster tools
#include "IClusterTool.hpp"
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "Cluster.h"

/**   
* @class FuzzyClusterTool
*
* Performs fuzzy clustering on all CAL hits
*
*
* $Header$
*/


class FuzzyClusterTool : public Cluster {

public:
    
    //! constructor
    FuzzyClusterTool( const std::string& type, const std::string& name, const IInterface* parent);
    //! destructor
     ~FuzzyClusterTool() {}; 
    
     StatusCode initialize();

        
/*!Performs the reconstruction, creates one or more CalCluster objects 
 * and stores there the following results: 
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

    //!  fuzzy cluster tool name
    std::string m_clusterToolName;

    //! pointer to the fuzzy cluster tool
    IClusterTool * m_fuzzyClusterTool;

    //! fuzzy analysis number
    int m_clusterSetNo;

};

#endif



