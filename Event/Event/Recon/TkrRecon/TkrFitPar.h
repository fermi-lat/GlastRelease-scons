//----------------------------------------
//
//      Part of the KALMAN Filter Objects Declarations
//
//      due to Brian Baughman 18/06/02
//      
//      
//----------------------------------------

#ifndef _TkrFitPar_H
#define _TkrFitPar_H 1

//#define TIMETRIAL 1
#undef TIMETRIAL

#include "CLHEP/Matrix/Vector.h"



/** 
 * @class TkrFitPar
 *
 * @brief Contains the Track Parameters
 *
 * @author Brian Baughman
 *
 * $Header$
 */

class Ray; 

namespace Event { // Namespace

  class TkrFitPar
  {
    // TkrFitPar Parameters arranged as a 4 - HepVector: 
    //       X intercept, X slope, Y intercept, Y slope
      
  private:
 
    double m_xpos;
    double m_xslope;
    double m_ypos;
    double m_yslope;

    friend class TkrFitMatrix;

  public:
    
    // Constructors: 
    // Creates a null parameter vector
    TkrFitPar ();

    TkrFitPar(const TkrFitPar &p);
      
    // Creat from a HepVector
    TkrFitPar(const HepVector &p);
#ifndef TIMETRIAL
	
    // Create from a 3D trajectory
    TkrFitPar (const Ray &); 
    TkrFitPar (const Ray *); 

#endif TIMETRIAL
    // Create with explicit values
    TkrFitPar(double ax, double sx, double ay,double sy);

       
    // Access methods for individual fit parameters
    double TkrFitPar::getXPosition() const;
    double TkrFitPar::getXSlope()    const;
	   
    double TkrFitPar::getYPosition() const;
    double TkrFitPar::getYSlope()    const;
	       
	       
    //operator overload for element access TkrFitPar(element)
    //!NOTE index starts with (1). In case of index error returns (1)
    const double & TkrFitPar::operator() (const int &element) const;
		 
		 
    //!NOTE index starts with (1). In case of index error returns (1)
    double & TkrFitPar::operator() (const int &element);
	    
    //operator overload for element access TkrFitPar(element)
    //!NOTE index starts with [0]. In case of index error returns [0]
    const double & TkrFitPar::operator[] (const int &element) const;
		 
		 
    //!NOTE index starts with [0]. In case of index error returns [0]
    double & TkrFitPar::operator[] (const int &element);


    //operator overload for TkrFitPar+TkrFitPar
    const TkrFitPar operator +(const TkrFitPar& A) const;
	      
    //operator overload for TkrFitPar-TkrFitPar
    const TkrFitPar operator -(const TkrFitPar& A) const;
		
    //operator overload for TkrFitPar*TkrFitPar
    const double operator *(const TkrFitPar& A) const;

  };
    


}; //Namespace

#endif _TkrFitPar_H

