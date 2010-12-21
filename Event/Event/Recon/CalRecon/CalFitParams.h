
#ifndef CalFitParams_H
#define CalFitParams_H


/** 
* @class CalFitParams
*
* @brief Gaudi TDS class to store the calorimeter fit parameters and associated
* covariance matrix.
*
* This is an extension of the CalParams base class meant to include all the
* quantities calculated in the fit to the calorimeter cluster cetroid/direction.
*
* Note that in the current implementation the m_energy and m_eneError class members
* are essentially unused, but it's not unconceivable that in the future the fit
* to the cluster might provide an energy measurement, so I decided not to refactor
* out these two quantities.
*
* Whenever this class is changed, the changes should be propagated to the
* related files on the ROOT side:
* - reconRootData/reconRootData/CalFitParams.h
* - reconRootData/src/CalFitParams.cxx
* 
* @author Bill Atwood, Tracy Usher, Luca Baldini
*
*/

#include <iostream>
#include "CalParams.h"


namespace Event { //Namespace Event


  class CalFitParams : public CalParams
  {
  public:
    /// Default constructor
    CalFitParams() { clear() ; }
 
    /// Direct construction from all the elements.
    CalFitParams(double energy, double eneError,
		 const Point& centroid, const CLHEP::HepMatrix& centroidErr,
		 const Vector& axis, const CLHEP::HepMatrix& axisErr,
		 int numFitLayers, double chiSquare);

    /// And even more parameters (reflecting the old-fashioned way CalParams constructor).
    CalFitParams(double energy, double eneError,
		 double xCntrd, double yCntrd, double zCntrd,
		 double cntdxx, double cntdxy, double cntdxz,
		 double cntdyy, double cntdyz, double cntdzz,
		 double xAxis,  double yAxis,  double zAxis,
		 double axsdxx, double axsdxy, double axsdxz,
		 double axsdyy, double axsdyz, double axsdzz,
		 int numFitLayers, double chiSquare);
 
    /// Convenience constructor (given Points, Vectors and matrices) for backward
    /// compatibility.
    /// The energy and energy errors are set to unphysical values. Also note that the
    /// order of the parameters refelects the old interface and is different with respect
    /// to the two previous constructors.
    CalFitParams(int numFitLayers, double chiSquare,
                 const Point& centroid, const CLHEP::HepMatrix& centroidErr,
                 const Vector& axis, const CLHEP::HepMatrix& axisErr);

    /// Convenience direct constructor from all the elements (the old fashioned way)
    /// (see the comment to the previous constructor).
    CalFitParams(int numFitLayers, double chiSquare,
                 double xCntrd, double yCntrd, double zCntrd,
                 double cntdxx, double cntdxy, double cntdxz,
		 double cntdyy, double cntdyz, double cntdzz,
                 double xAxis,  double yAxis,  double zAxis,
                 double axsdxx, double axsdxy, double axsdxz,
		 double axsdyy, double axsdyz, double axsdzz);

    /// Destructor;
    ~CalFitParams() {}
    
    /// reset method
    void clear() ;
    /// Part of the reset code specific to the CalFitParams class
    /// (as opposed to the base class CalParams).
    void clearFitParams();
    
    /// Retrieve parameters...
    inline int getFitLayers()            const {return m_nFitLayers;}
    inline double getChiSquare()         const {return m_chiSquare;}

    /// Set parameters.
    inline void setFitLayers(int numFitLayers) {m_nFitLayers = numFitLayers;}
    inline void setChiSquare(double chiSquare) {m_chiSquare  = chiSquare;}

    std::ostream& fillStream( std::ostream& s ) const;
    friend std::ostream& operator<< ( std::ostream& s, const CalFitParams& obj ) 
    {
      return obj.fillStream(s);
    }

private:

    /// Parameters associated with the fit
    int    m_nFitLayers;
    double m_chiSquare;
  };


}; //Namespace Event


#endif

