
#ifndef __IEnergyCorr_H
#define __IEnergyCorr_H 1

#include "GaudiKernel/IAlgTool.h"
#include "Event/Recon/CalRecon/CalCluster.h"
#include "geometry/Vector.h"

/**   
* @class IEnergyCorr
*
* base 
*
*
* Performs high level energy corrections
*
* The reconstruction here uses CalXtalRecCol to produce a CalEnergyCorrCol.
*  It evaluates the barycenter for each layer using the coordinates stored in
*  the CalXtalRecCol,
*
* $Header$
*/

static const InterfaceID IID_IEnergyCorr("IEnergyCorr", 1 , 0);

class IEnergyCorr : virtual public IAlgTool {

public:
	IEnergyCorr() {;};
    // retrieve interface ID
    static const InterfaceID& interfaceID() { return IID_IEnergyCorr; }
    //! constructor
    //! destructor
    virtual ~IEnergyCorr() {}; 
    
    virtual StatusCode initialize()=0;

        
/*!Performs the reconstruction, creates one CalEnergyCorr object and stores
 * there the following results: 
 * - Energy per layer is computed and stored in CalEnergyCorr in MeV
 * - Barycenter per layer is also computed and stored in CalEnergyCorr
 */        

    virtual StatusCode doEnergyCorr(double eTotal, Event::CalCluster* cluster)=0;
    virtual StatusCode execute()=0;

    virtual StatusCode finalize()=0; 

};

#endif



