
#ifndef __IEnergyCorr_H
#define __IEnergyCorr_H 1

#include "GaudiKernel/IAlgTool.h"
#include "Event/Recon/CalRecon/CalCluster.h"
#include "geometry/Vector.h"
#include "CalClusteringData.h"

/**   
* @class IEnergyCorr
*
* base class for energy leakage corrections 
*
*
* $Header$
*/

static const InterfaceID IID_IEnergyCorr("IEnergyCorr", 1 , 0) ;

class IEnergyCorr : virtual public IAlgTool {

public:

    // retrieve interface ID
    static const InterfaceID & interfaceID() { return IID_IEnergyCorr ; }
    
    //! constructor
	IEnergyCorr() {}
    //! destructor
    virtual ~IEnergyCorr() {}
    
    //! Worker function for calculating corrections
	/*! Performs the reconstruction, creates one CalEnergyCorr object and stores
     *  there the following results: 
     *  - Energy per layer is computed and stored in CalEnergyCorr in MeV
     *  - Barycenter per layer is also computed and stored in CalEnergyCorr
     */        
    virtual StatusCode doEnergyCorr( const CalClusteringData *, Event::CalCluster * ) =0 ;
    
    virtual double getEnergyCorr() =0 ;

} ;

#endif



