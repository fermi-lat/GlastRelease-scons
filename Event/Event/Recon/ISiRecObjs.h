
#ifndef __ISIRECOBJS_H
#define __ISIRECOBJS_H 1

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/DataObject.h"
#include "geometry/Point.h"
#include "geometry/Vector.h"

//----------------------------------------------
//
//   ISiRecObjs
//
//   Interface for Transient Storage Data
//----------------------------------------------
//   Access for very basic aspects of SiRecObjs 
//   mainly for DOCA and display
//   
//----------------------------------------------
//             Ian Gable, Victoria 06/21/2001
//----------------------------------------------

extern const CLID& CLID_SiRecObjs;
//const static CLID CLID_SiRecObjs = 254;

/*!
ISiRecObjs striped down interface for the SiRecObjs class
*/
//##########################################################
class ISiRecObjs : public DataObject
//##########################################################
{

public:

	//! constructor - ini the list
	ISiRecObjs()  {;};
	//! destructor - delete the pointers to the data
	virtual ~ISiRecObjs() {}
	
	//! returns number of gammas
	virtual int numGammas() const=0;
	//! returns number of GFparticles
	virtual int numParticles() const=0;

    //! Get the X slope of the ith GFparticle
    virtual double getXGFparticleSlope(int i) =0;
    //! Get the Y slope of the ith GFparticle
	virtual double getYGFparticleSlope(int i) =0;
    //! get the vextex of the ith Gamma
    virtual Point getGammaVertex(int i) =0;
    //! get the direction vector of the ith Gamma
    virtual Vector getGammaDirection(int i) =0;
    //! Get the X slope of the ith GFgamma
    virtual double getXGFgammaSlope(int i) =0;
    //! Get the Y slope of the ith GFgamma
    virtual double getYGFgammaSlope(int i)  =0;

private:


};
      
#endif
