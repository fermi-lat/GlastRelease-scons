#ifndef ICsIClusterList_H
#define ICsIClusterList_H

#include <vector>
#include "geometry/Point.h"
#include "geometry/Vector.h"
#include "GaudiKernel/DataObject.h"

extern const CLID& CLID_CalClusterList;


//----------------------------------------------
//
//   ICsICluster
//
//   Interface for Transient Storage Data
//----------------------------------------------
//   Invented for the use of TkrRecon
//----------------------------------------------
//            LSR  09/14/2001
//----------------------------------------------


//#####################################
class ICsICluster
//#####################################
{
protected:
	//! constructor
	ICsICluster() {};

public:
	//! Destructor
    virtual ~ICsICluster() {};

    // TkrRecon uses only Energy and position so far, but might as
    //   well keep all the access functions for now.

	// access
	virtual double energySum()        const = 0;
	virtual double energyLeak()       const = 0;
	virtual double energyCorrected()  const = 0;
	virtual double getEneLayer(int i) const = 0;
	virtual const  Vector& getPosLayer(int i)         const = 0;
	virtual const  std::vector<double>& getEneLayer() const = 0;
	virtual const  std::vector<Vector>& getPosLayer() const = 0;
	virtual const  std::vector<Vector>& getRmsLayer() const = 0;
	virtual double getRmsLong()		  const = 0;
	virtual double getRmsTrans()	  const = 0;
    virtual double getTransvOffset()  const = 0;

	virtual Point  position()         const = 0;
	virtual Vector direction()        const = 0;
	virtual double getFitEnergy()	  const = 0;
	virtual double getProfChisq()	  const = 0;
	virtual double getCsiAlpha()	  const = 0;
	virtual double getCsiLambda()	  const = 0;
	virtual double getCsiStart()	  const = 0;
	// operations
	virtual void   writeOut() const = 0;

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
class ICsIClusterList : public DataObject
//#####################################
{
protected:
    ICsIClusterList() {};

public:
    virtual ~ICsIClusterList() {};

	// GAUDI members to be use by the converters
	static const CLID& classID() {return CLID_CalClusterList;}
	virtual const CLID& clID() const {return classID();}

	// access
	virtual int num()                   const = 0;
	virtual ICsICluster* Cluster(int i) const = 0;

};
#endif	






