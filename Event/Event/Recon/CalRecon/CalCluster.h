#ifndef CalCluster_H
#define CalCluster_H

#include <vector>
#include "geometry/Point.h"
#include "geometry/Vector.h"
#include "GaudiKernel/DataObject.h"
#include "GaudiKernel/MsgStream.h"

// #include "CalRecon/CalDisplay.h"

extern const CLID& CLID_CalClusterCol;


//----------------------------------------------
//
//   CalCluster
//
//   Transient Storage Data
//----------------------------------------------
//   It contains the high level data for the Calorimeter
//----------------------------------------------
//             J.A Hernando, Santa Cruz 02/29/00
//----------------------------------------------

//! High level data for the calorimeter
/*! Transient storage of the results of the reconstruction performed
 *  in CalClustersAlg. It contains the data from one cluster in the
 *  calorimeter. 
 *  
 *  \author Alexandre Chehtman
 *  \author Regis Terrier
 *  \author Jose Angel Hernando
 *
 * \b Revisions:
 *  - 10/17/00    RT   Comments added
 *  - 08/20/00    RT   m_leakEnergy added
 *  - 02/29/00    JAH  first implementation
 */

namespace Event 

{

//#####################################
class CalCluster 
//#####################################
{
public:

	//! constructor
	CalCluster(double e, Point p);

	//! Destructor
	~CalCluster() {}
	void initialize(double eleak, std::vector<double> eneLayer, std::vector<Vector> pLayer,
					std::vector<Vector> rmsLayer, double rms_long, double rms_trans,
		            Vector caldir, double calTransvOffset)
	{
	m_eneLayer = eneLayer;
	m_leakEnergy = eleak;
	m_pLayer = pLayer;
	m_rmsLayer = rmsLayer;	
	m_rmslong = rms_long;
	m_rmstrans= rms_trans;
	m_direction = caldir;
	m_transvOffset = calTransvOffset;
	}

    void initProfile(double fit_energy, double ki2, double fit_start,
		double fit_alpha, double fit_lambda)
	{
	m_fitEnergy = fit_energy;
	m_ProfChisq = ki2;
	m_CsiAlpha = fit_alpha;
	m_CsiLambda = fit_lambda;
	m_start  = fit_start;
	}


	// access
	double getEnergySum()        const {return m_energySum;}
	double getEnergyLeak()       const {return m_leakEnergy;}
	double getEnergyCorrected()  const {return m_energyCorrected;}
	double getEneLayer(int i) const {return m_eneLayer[i];}
	const Vector& getPosLayer(int i) const {return m_pLayer[i];}
	const std::vector<double>& getEneLayer() const {return m_eneLayer;}
	const std::vector<Vector>& getPosLayer() const {return m_pLayer;}
	const std::vector<Vector>& getRmsLayer() const {return m_rmsLayer;}
	double getRmsLong()		  const {return m_rmslong;}
	double getRmsTrans()	  const {return m_rmstrans;}
    double getTransvOffset()  const {return m_transvOffset;}

	Point getPosition()          const {return m_position;}
	Vector getDirection()        const {return m_direction;}
	double getFitEnergy()	  const {return m_fitEnergy;}
	double getProfChisq()	  const {return m_ProfChisq;}
	double getCsiAlpha()	  const {return m_CsiAlpha;}
	double getCsiLambda()	  const {return m_CsiLambda;}
	double getCsiStart()	  const {return m_start;}
	// operations
	void writeOut(MsgStream& stream) const;

protected:

	virtual void ini();

private:

	//! Total measured energy in the calorimeter
	double m_energySum;
	//! Leakage corrected energy using correlation method ( for E> several GeV)
	double m_leakEnergy;
	//! corrected energy not used ( yet )
	double m_energyCorrected;
	//! Energy per layer in MeV
	std::vector<double> m_eneLayer;
	//! Barycenter position in each layer
	std::vector<Vector> m_pLayer;
	//! RMS of energy deposition in each layer
	std::vector<Vector> m_rmsLayer;
	//! RMS of longitudinal position measurement
	double m_rmslong;
	//! RMS of transverse position measurement
	double m_rmstrans;
    //! Transvers offset of calorimeter position measurement
    double m_transvOffset;
	
	//! fitted energy ( for E>10 GeV)
	double m_fitEnergy;
	//! Chisquare of the fit ( not a real Chisquare)
	double m_ProfChisq;
	//! Alpha parameter used in the fit
	double m_CsiAlpha;
	//! Lambda parameter used in the fit 
	double m_CsiLambda;
	//! Fitted starting point of the shower (physical meaning is not clear)
	double m_start;

	Point m_position;
	Vector m_direction;
};


//! High level data for the calorimeter
/*! Transient storage of the results of the reconstruction performed
 *  in CalClustersAlg. It contains the the list ofall the clusters in the
 *  calorimeter. 
 *  
 * \warning there is no clustering in the calorimeter up to now. There is
 *  just one Cluster stored containing the information for the whole event
 *
 *  \author Jose Angel Hernando
 *
 * \b Revisions:
 * - 02/29/00     JAH    first implementation
 */

//#####################################
class CalClusterCol : public DataObject, public std::vector<CalCluster*>
//#####################################
{
public:

	CalClusterCol() { clear();}
	~CalClusterCol() { delClusters();}


	// GAUDI members to be use by the converters
	static const CLID& classID() {return CLID_CalClusterCol;}
	virtual const CLID& clID() const {return classID();}
	
	
	void add(CalCluster* cl) {push_back(cl);}

	// access
	int num()                  const {return size();}
	CalCluster* getCluster(int i) const {return operator[](i);}

	//operations
	void delClusters();

	virtual void writeOut(MsgStream& stream) const;

protected:

	virtual void ini();

};

}

#endif	






