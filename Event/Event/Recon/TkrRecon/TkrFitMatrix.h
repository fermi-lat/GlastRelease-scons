//----------------------------------------
//
//      Kalman Filter Objects Declarations
//
//      Original due to Brian Baughman 17-6-2002
//      rewritten from Jose's code
//      
//----------------------------------------

#ifndef _TkrFitMatrix_H
#define _TkrFitMatrix_H 1

//#define TIMETRIAL 1
#undef TIMETRIAL

#ifdef TIMETRIAL
#include "TkrFitPar.h"

#else /* Standard Gleam */
#include "Event/Recon/TkrRecon/TkrFitPar.h"

#endif /* TIMETRIAL */

#include "CLHEP/Matrix/Matrix.h"
#include <cfloat>




/** 
 * @class TkrFitHit
 *
 * @brief Contains the Track Covariance matrix
 *
 * @author Brian Baughman
 *
 * $Header$
 */

namespace Event { //Namespace
  
  class TkrFitMatrix 
  {
    // Kalman matrices: block diagonal 2 - 2x2's as a 4x4 matrix. 
    // Upper block: X - projection; Lower Block: Y - projection
    // Note:  HepMatrix class indexs from 1 not 0
    
  private:

    double m_11;
    double m_12;
    double m_13;     
    double m_14;

    double m_21;
    double m_22;
    double m_23;
    double m_24;
 
    double m_31;
    double m_32;
    double m_33;
    double m_34;
	 
    double m_41;
    double m_42;
    double m_43;
    double m_44;
		 
    double m_zero;


		     
  public:
    // Constructors:
    // A null matrix
    TkrFitMatrix();
    // A matrix
    TkrFitMatrix(
		 double e_11, double e_12, double e_13, double e_14,
		 double e_21, double e_22, double e_23, double e_24,
		 double e_31, double e_32, double e_33, double e_34,
		 double e_41, double e_42, double e_43, double e_44
		 );
    // constructor from CLHEP Matrix
    TkrFitMatrix(const HepMatrix &A);
    
    // constructor from TkrFitMatrix
    TkrFitMatrix(const TkrFitMatrix &A);

    // Special constructor to produce propagation matrix F
    TkrFitMatrix(double step);

    // Special constructor to  produce matrix H
    TkrFitMatrix(int one);
		   
    // Access to elements of the covariance matrix
    inline double TkrFitMatrix::getcovX0X0() const {return m_11;}
    inline double TkrFitMatrix::getcovSxSx() const {return m_22;}
    inline double TkrFitMatrix::getcovX0Sx() const {return m_12;}
    inline double TkrFitMatrix::getcovSxX0() const {return m_21;}
			   
    inline double TkrFitMatrix::getcovY0Y0() const {return m_33;}
    inline double TkrFitMatrix::getcovSySy() const {return m_44;}
    inline double TkrFitMatrix::getcovY0Sy() const {return m_34;}
    inline double TkrFitMatrix::getcovSyY0() const {return m_43;}
				   
    inline double TkrFitMatrix::getcovX0Y0() const {return m_13;}
    inline double TkrFitMatrix::getcovX0Sy() const {return m_14;}
    inline double TkrFitMatrix::getcovSyX0() const {return m_41;}
						
    inline double TkrFitMatrix::getcovY0X0() const {return m_31;}
    inline double TkrFitMatrix::getcovY0Sx() const {return m_32;}
    inline double TkrFitMatrix::getcovSxY0() const {return m_23;}
					     
    inline double TkrFitMatrix::getcovSxSy() const {return m_24;}
    inline double TkrFitMatrix::getcovSySx() const {return m_42;}
						 

    // function for calculating the transpose of the TkrFitMatrix in question
    TkrFitMatrix TkrFitMatrix::T() const;
		 
    // function for calculating the inverse of the TkrFitMatrix in question
    // NOTE: ierr = 0 is clean run else is NOT inverted
    TkrFitMatrix TkrFitMatrix::inverse(int &ierr) const;
		   
    // function for calculating the inverse of TkrFitMatrix REDEFINES object
    // NOTE: ierr = 0 is clean run else is NOT inverted
    void TkrFitMatrix::invert(int &ierr);

    //operator overload for access/definition TkrFitMatrix(row,column)
    //!NOTE index starts with (1,1)
    const double & TkrFitMatrix::operator() 
      (const int &row, const int &column) const;

    double & TkrFitMatrix::operator() (const int &row, const int &column);
    
    //operator overload for TkrFitMatrix+TkrFitMatrix
    const TkrFitMatrix operator +(const TkrFitMatrix& A) const;
						   
    //operator overload for TkrFitMatrix-TkrFitMatrix
    const TkrFitMatrix operator -(const TkrFitMatrix& A) const;
						     

    //operator overload for TkrFitMatrix*TkrFitMatrix
    //!NOTE: Assumed LH multiply
    const TkrFitMatrix operator *(const TkrFitMatrix& A) const;
	
    //operator overload for TkrFitMatrix*TkrFitPar
    const TkrFitPar operator *(const TkrFitPar& A) const;
    
    //operator overload for TkrFitMatrix*double
    const TkrFitMatrix operator *(const double& A) const;	   
  };

}; //Namespace

#endif _TkrFitMatrix_H






