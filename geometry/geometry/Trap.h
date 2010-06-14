// $Id$
//


#ifndef __TRAP_H
#define __TRAP_H

#include "geometry/Hexahedron.h"
 
/**
	Trap is the same as the GEANT volume of the same name, but note that
	GEANT uses half-lengths as parameters, whereas we use full
	lengths.<BR>
	A Trap is defined by the  6 planes which bound it.
	It is a generalized trapezoidal volume.
	The GEANT specification for this shape takes 11 parameters as
	described below, but, as noted by the authors, only 9 of them are 
	independent.
	Here, we use only the 9 independent parameters, eliminating the ones
	indicated below.  <BR>
	The 11 GEANT parameters are
        <OL>
	 <LI> the full length in z
	 <LI> the polar angle of the center of the -z face to the center of 
		the +z face
	 <LI> the azimuthal angle of the center of the -z face to the center
		of the +z face
		Since this parameter is not independent of overall rotations,
		THIS PARAMETER IS not used!
	 <LI> the full length in y of the -z face
	 <LI> the full length in x of the low y edge of the -z face
	 <LI> the full length in x of the high y edge of the -z face
	 <LI> the angle with respect to the y-axis from the center of the
		low y edge to the center of the high y edge of the -z face
	 <LI> the full length in y of the +z face
	<LI> the full length in x of the low y edge of the +z face
	<LI> the full length in x of the high y edge of the +z face
	<LI> the angle with respect to the y-axis from the center of the
		low y edge to the center of the high y edge of the +z face
		Since this parameter must be equal to 7) in order for all
		faces to form planes,
		THIS PARAMETER IS not used!
       </OL>

*/
class Trap: public Hexahedron
{

public:
	Trap( double h, double th, 
	      double w1, double bl1, double tl1, double th1,
	      double w2, double bl2, double tl2 ); 


	double height() const ;
	double theta() const;
	double width1() const; 
	double lengthLow1() const; 
	double lengthHigh1() const ;
	double theta1() const ; 
	double width2() const ;
	double lengthLow2() const ; 
	double lengthHigh2() const ;
	// functions that should return the original parameters   

	virtual const char *nameOf() const { return "Trap"; }
	virtual void printOn( std::ostream& os = std::cout ) const;

};

#endif

