#ifndef CsIClusterList_H
#define CsIClusterList_H

#include <vector>
#include "geometry/Point.h"
#include "geometry/Vector.h"
// #include "Event/trsDataVI.h"
#include "GaudiKernel/DataObject.h"

#include "CalRecon/CalDisplay.h"

extern const CLID& CLID_CalClusterList;


//----------------------------------------------
//
//   CsICluster
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


//#####################################
class CsICluster
//#####################################
{
public:

	//! constructor
	CsICluster(double e, Point p);

	//! Destructor
	~CsICluster() {}

	void setDirection(Vector v)   {m_direction = v;}

	//! Set energy corrected
	/*! not used for the moment
	 *  Energy sum is stored */
	void setEnergyCorrected(double e) {m_energyCorrected = e;}
	//! Set energy per layer
	void setEneLayer(std::vector<double> v){m_eneLayer = v;}
	//! Set barycenter position for each layer
	void setPosLayer(std::vector<Vector> v){m_pLayer = v;}
	//! Set rms of energy deposition for each layer
	void setRmsLayer(std::vector<Vector> v){m_rmsLayer = v;}
	//! Set Longitudinal RMS
	void setRmsLong(double r) {m_rmslong=r;}
	//! Set transverse RMS
	void setRmsTrans(double r) {m_rmstrans=r;}
	//! Set energy corrected via CalClustersAlg::Leak()
	void setEneLeak(double e) {m_leakEnergy = e;}
	//! Set fitted energy form CalClustersAlg::Profile()
	void setFitEnergy(double e) { m_fitEnergy = e;}
	//! Set chi square of profile fitting
	void setProfChisq(double k) { m_ProfChisq = k;}
	//! Set alpha parameter used in the fit
	void setCsiAlpha(double a) { m_CsiAlpha =a;}
	//! Set lambda parameter used in the fit
	void setCsiLambda(double l) { m_CsiLambda = l;}
	//! Set the fitted starting point
	void setCsiStart(double s) { m_start = s;}
    //! Set the transverse offset of calorimeter position measurement
    void setTransvOffset (double offset) {m_transvOffset = offset;}

	// access
	double energySum()        const {return m_energySum;}
	double energyLeak()       const {return m_leakEnergy;}
	double energyCorrected()  const {return m_energyCorrected;}
	double getEneLayer(int i) const {return m_eneLayer[i];}
	const Vector& getPosLayer(int i) const {return m_pLayer[i];}
	const std::vector<double>& getEneLayer() const {return m_eneLayer;}
	const std::vector<Vector>& getPosLayer() const {return m_pLayer;}
	const std::vector<Vector>& getRmsLayer() const {return m_rmsLayer;}
	double getRmsLong()		  const {return m_rmslong;}
	double getRmsTrans()	  const {return m_rmstrans;}
    double getTransvOffset()  const {return m_transvOffset;}

	Point position()          const {return m_position;}
	Vector direction()        const {return m_direction;}
	double getFitEnergy()	  const {return m_fitEnergy;}
	double getProfChisq()	  const {return m_ProfChisq;}
	double getCsiAlpha()	  const {return m_CsiAlpha;}
	double getCsiLambda()	  const {return m_CsiLambda;}
	double getCsiStart()	  const {return m_start;}
	// operations
	void writeOut() const;

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
class CsIClusterList : public DataObject
//#####################################
{
public:

	CsIClusterList() { m_calDisp = 0;clear();}
	~CsIClusterList() { if(m_calDisp) m_calDisp->clearClusterDisp();clear();}


	// GAUDI members to be use by the converters
	static const CLID& classID() {return CLID_CalClusterList;}
	virtual const CLID& clID() const {return classID();}
	
	
	void add(CsICluster* cl) {m_CsIClustersList.push_back(cl);}

	// access
	int num()                  const {return m_CsIClustersList.size();}
	CsICluster* Cluster(int i) const {return m_CsIClustersList[i];}

	//operations
	virtual void clear();
	virtual void make() {}

	virtual void writeOut() const;
    void setCalDisplay(CalDisplay* calDisp) {m_calDisp = calDisp;}

protected:

	virtual void ini();

private:

	std::vector<CsICluster*> m_CsIClustersList;

    CalDisplay* m_calDisp;
};
#endif	






