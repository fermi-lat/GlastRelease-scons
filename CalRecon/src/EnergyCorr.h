
#ifndef __EnergyCorr_H
#define __EnergyCorr_H 1

#include "IEnergyCorr.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/ISvcLocator.h"
#include "Event/Recon/CalRecon/CalCluster.h"

/**   
* @class EnergyCorr
*
* Algorithm for reconstruction of energy and direction of incident particle
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


class EnergyCorr : virtual public IEnergyCorr {
	
public:
    
    //! constructor
	
    EnergyCorr() {};
    //! destructor
    virtual ~EnergyCorr() {}; 
    
    virtual StatusCode initialize()  {
		return StatusCode::SUCCESS;
	};
	
	
	/*!Performs the reconstruction, creates one CalEnergyCorr object and stores
	* there the following results: 
	* - Energy per layer is computed and stored in CalEnergyCorr in MeV
	* - Barycenter per layer is also computed and stored in CalEnergyCorr
	*/        
    virtual StatusCode doEnergyCorr(double eTotal, Event::CalCluster* cluster) 
	{return StatusCode::SUCCESS;};
	
    virtual StatusCode execute() {return StatusCode::SUCCESS;};
	
    virtual StatusCode finalize() {return StatusCode::SUCCESS;};
    
    virtual double getEnergyCorr() {return m_energyCorr;};

    virtual double getTrackSlope() {return m_slope;};

    virtual void setTrackSlope(double slope) {m_slope=slope;};

    virtual int getNLayers() {return m_nLayers;};

protected:

    virtual void setEnergyCorr(double energy) { m_energyCorr = energy;};
    virtual void setNLayers(int n) { m_nLayers = n;};
	
private:
	
    // energy correction
    double m_energyCorr;
    // number of layers in CAL
    int m_nLayers;
    // slope of track
    double m_slope;
};
#endif
	
	
	
