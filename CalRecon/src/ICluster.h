
#ifndef __ICluster_H
#define __ICluster_H 1

#include "GaudiKernel/IAlgTool.h"
#include "Event/Recon/CalRecon/CalXtalRecData.h"

/**   
* @class ICluster
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

static const InterfaceID IID_ICluster("ICluster", 1 , 0);

class ICluster : virtual public IAlgTool {

public:
	ICluster() {;};
    // retrieve interface ID
    static const InterfaceID& interfaceID() { return IID_ICluster; }
    //! constructor
    //! destructor
    virtual ~ICluster() {}; 
    
    virtual StatusCode initialize()=0;

        
/*!Performs the reconstruction, creates one CalCluster object and stores
 * there the following results: 
 * - Energy per layer is computed and stored in CalCluster in MeV
 * - Barycenter per layer is also computed and stored in CalCluster
 */        

    virtual StatusCode findClusters(Event::CalXtalRecCol* calXtalRecCol)=0;
    
	virtual StatusCode execute()=0;

    virtual StatusCode finalize()=0;
    

    
    
    
};

#endif



