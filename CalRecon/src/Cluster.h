
#ifndef __Cluster_H
#define __Cluster_H 1

#include "ICluster.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/ISvcLocator.h"
#include "geometry/Vector.h"
#include "Event/Recon/CalRecon/CalXtalRecData.h"
#include "Event/Recon/CalRecon/CalCluster.h"

/**   
* @class Cluster
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


class Cluster : virtual public ICluster {
	
public:
    
    //! constructor
	
    Cluster() {};
    //! destructor
    virtual ~Cluster() {}; 
    
    virtual StatusCode initialize()  {
		return StatusCode::SUCCESS;
	};
	
	
	/*!Performs the reconstruction, creates one CalCluster object and stores
	* there the following results: 
	* - Energy per layer is computed and stored in CalCluster in MeV
	* - Barycenter per layer is also computed and stored in CalCluster
	*/        
    virtual StatusCode findClusters(Event::CalXtalRecCol* calXtalRecCol) 
	{return StatusCode::SUCCESS;};
	
    virtual StatusCode execute() {return StatusCode::SUCCESS;};
	
    virtual StatusCode finalize() {return StatusCode::SUCCESS;};
    
	virtual Event::CalXtalRecCol* getRecCol() {return m_calXtalRecCol;};

	virtual Event::CalClusterCol* getClusterCol() {return m_calClusterCol;};
	
	virtual double getEnergy() {return m_ene;};
    
	virtual Vector getClusterPosition() {return m_pCluster;};
    
	virtual std::vector<double> getEnergyPerLayer() 
	{return m_eneLayer;};
	
	virtual std::vector<Vector> getPositionPerLayer() 
	{return m_pLayer;};
	
	virtual std::vector<Vector> getRmsPerLayer() 
	{return m_rmsLayer;};
	
	virtual double getRmsTrans() {return m_rms_trans;};
	
	virtual double getRmsLong() {return m_rms_long;};

	virtual void setClusterCol(Event::CalClusterCol* calClusterCol)
			{m_calClusterCol = calClusterCol;};

	
	virtual void setData(Event::CalXtalRecCol* calXtalRecCol,
		double energy,
		std::vector<double> EnergyVector,
		Vector pos,
		std::vector<Vector> PositionVector,
		std::vector<Vector> RmsVector,
		double rmsTrans,
		double rmsLong) {
		setRecCol(calXtalRecCol);
		setEnergy(energy);
		setClusterPosition(pos);
		setEnergyPerLayer(EnergyVector);
		setPositionPerLayer(PositionVector);
		setRmsPerLayer(RmsVector);
		setRmsTrans(rmsTrans);
		setRmsLong(rmsLong);	
	}
	
protected:
	
	virtual void setRecCol(Event::CalXtalRecCol* calXtalRecCol)
	{m_calXtalRecCol = calXtalRecCol;};
	virtual void setClusterPosition(Vector pos) {m_pCluster = pos;};
	virtual void setEnergy(double energy) {m_ene = energy;};
	virtual void setEnergyPerLayer(std::vector<double> EnergyVector)
	{m_eneLayer = EnergyVector;};
	virtual void setPositionPerLayer(std::vector<Vector> PositionVector)
	{m_pLayer = PositionVector;};
	virtual void setRmsPerLayer(std::vector<Vector> RmsVector)
	{m_rmsLayer = RmsVector;};
	virtual void setRmsTrans(double rmsTrans) {m_rms_trans = rmsTrans;};
	virtual void setRmsLong(double rmsLong) {m_rms_long = rmsLong;};
	
	
private:
	
	
	//! reconstructed data for crystals, the input of the reconstruction
	Event::CalXtalRecCol* m_calXtalRecCol;
	
	// total energy in calorimeter
	double m_ene;                 
	
	// cluster position
	Vector m_pCluster;            
	
	// energy per layer
	std::vector<double> m_eneLayer;   
	
	// Vector of average position per layer
	std::vector<Vector> m_pLayer;    
	
	// Vector of quadratic spread per layer
	std::vector<Vector> m_rmsLayer;  
	
	// transverse rms
	double m_rms_trans;
	// longitudinal rms
	double m_rms_long;
	
	Event::CalClusterCol* m_calClusterCol;
	};


#endif
	
	
	
