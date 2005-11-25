
#ifndef CalParams_H
#define CalParams_H

#include "geometry/Point.h"
#include "geometry/Vector.h"
#include "CLHEP/Matrix/Matrix.h"

/** 
* @class CalParams
*
* @brief Gaudi TDS class to store the calorimeter parameters and associated covariance
*        matrix. 
* 
* @author Bill Atwood, Tracy Usher
*
* $Header$
*/

#include <iostream>

namespace Event { //Namespace Event


class CalParams
{
public:
    /// Default constructor
    CalParams() { clear() ; }

    /// reset method
    void clear() ;
  
    /// Constructor given Points, Vectors and matrices
    CalParams(double energy, double eneError,
              const Point&  centroid, const HepMatrix& centroidErr,
              const Vector& axis,     const HepMatrix& axisErr);

    /// Direct construction from all the elements (the old fashioned way)
    CalParams(double energy, double eneError,
              double xCntrd, double yCntrd, double zCntrd,
              double cntdxx, double cntdxy, double cntdxz, double cntdyy, double cntdyz, double cntdzz,
              double xAxis,  double yAxis,  double zAxis,
              double axsdxx, double axsdxy, double axsdxz, double axsdyy, double axsdyz, double axsdzz);

    /// Copy constructor
    //CalParams (const CalParams& right);

   ~CalParams() {}

    /// Retrieve the energy
    inline double        getEnergy()       const {return m_energy;}
    inline double        getEnergyErr()    const {return m_eneError;}

    /// Retrieve the centroid position
    inline const Point&  getCentroid()     const {return m_clusterCentroid;}

    /// Errors in a HepMatrix
    HepMatrix            getCentroidErrs() const ;

    /// Direct access to errors
    inline double        getxPosxPos()     const {return m_cenxx; }
    inline double        getxPosyPos()     const {return m_cenxy; }
    inline double        getxPoszPos()     const {return m_cenxz; }
    inline double        getyPosyPos()     const {return m_cenyy; }
    inline double        getyPoszPos()     const {return m_cenyz; }
    inline double        getzPoszPos()     const {return m_cenzz; }

    /// Retrieve the centroid position
    inline const Vector& getAxis()         const {return m_clusterAxis;}

    /// Errors in a HepMatrix
    HepMatrix            getAxisErrs()     const ;

    /// Direct access to errors
    inline double        getxDirxDir()     const {return m_axisxx; }
    inline double        getxDiryDir()     const {return m_axisxy; }
    inline double        getxDirzDir()     const {return m_axisxz; }
    inline double        getyDiryDir()     const {return m_axisyy; }
    inline double        getyDirzDir()     const {return m_axisyz; }
    inline double        getzDirzDir()     const {return m_axiszz; }

    /// Define an ( ) operator (allows read/write - indexing from 1!!)
    //double& operator()(const int &i);
    //const double& operator()(const int &i) const;
    //double& operator()(const int &i, const int &j);
    //const double& operator()(const int &i, const int &j) const;

    /// Set the energy
    inline void setEnergy( double energy)       {m_energy   = energy;}
    inline void setEnergyErr( double energyErr) {m_eneError = energyErr;}

    /// Set parameters for centroid
    inline void setCentroid(const Point& pos)   {m_clusterCentroid = pos; }

    void        setCentroidErrs(const HepMatrix& errs);

    /// Set errors
    inline void setxDirxDir(double val)         {m_cenxx = val; }
    inline void setxDiryDir(double val)         {m_cenxy = val; }
    inline void setxDirzDir(double val)         {m_cenxz = val; }
    inline void setyDiryDir(double val)         {m_cenyy = val; }
    inline void setyDirzDir(double val)         {m_cenyz = val; }
    inline void setzDirzDir(double val)         {m_cenzz = val; }

    /// Set parameters for axis
    inline void setAxis(const Vector& axis)     {m_clusterAxis = axis; }

    void        setAxisErrs(const HepMatrix& errs);

    /// Set errors
    inline void setxPosxPos(double val)         {m_axisxx = val; }
    inline void setxPosyPos(double val)         {m_axisxy = val; }
    inline void setxPoszPos(double val)         {m_axisxz = val; }
    inline void setyPosySlp(double val)         {m_axisyy = val; }
    inline void setyPoszPos(double val)         {m_axisyz = val; }
    inline void setzPoszPos(double val)         {m_axiszz = val; }

    std::ostream& fillStream( std::ostream& s ) const;
    friend std::ostream& operator<< ( std::ostream& s, const CalParams& obj ) 
    {
	    return obj.fillStream(s);
    }

private:

    /// Energy and associated error
    double m_energy;
    double m_eneError;

    /// Cluster centroid
    Point  m_clusterCentroid;

    /// Error matrix for centroid stored in upper diagonal form
    double m_cenxx;     // Cov(1,1) = dx * dx 
    double m_cenxy;     // Cov(1,2) = Cov(2,1) = dx * dy
    double m_cenxz;     // Cov(1,3) = Cov(3,1) = dx * dz
    double m_cenyy;     // Cov(2,2) = dy * dy
    double m_cenyz;     // Cov(2,3) = Cov (3,2)= dy * dz
    double m_cenzz;     // Cov(3,3) = dz * dz

    /// Cluster axis, assumed to pass through centroid
    Vector m_clusterAxis;

    /// Error matrix for centroid stored in upper diagonal form
    double m_axisxx;    // Cov(1,1) = dx * dx 
    double m_axisxy;    // Cov(1,2) = Cov(2,1) = dx * dy
    double m_axisxz;    // Cov(1,3) = Cov(3,1) = dx * dz
    double m_axisyy;    // Cov(2,2) = dy * dy
    double m_axisyz;    // Cov(2,3) = Cov (3,2)= dy * dz
    double m_axiszz;    // Cov(3,3) = dz * dz
};


}; //Namespace Event


#endif

