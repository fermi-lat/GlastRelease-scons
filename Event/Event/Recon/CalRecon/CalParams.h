
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
    CalParams() {initDataMembers();}

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
    inline const double getEnergy()       const {return m_energy;}
    inline const double getEnergyErr()    const {return m_eneError;}

    /// Retrieve the centroid position
    inline const Point  getCentroid()     const {return m_clusterCentroid;}

    /// Errors in a HepMatrix
    const HepMatrix     getCentroidErrs() const;

    /// Direct access to errors
    inline const double getxPosxPos()     const {return m_cenxx; }
    inline const double getxPosyPos()     const {return m_cenxy; }
    inline const double getxPoszPos()     const {return m_cenxz; }
    inline const double getyPosyPos()     const {return m_cenyy; }
    inline const double getyPoszPos()     const {return m_cenyz; }
    inline const double getzPoszPos()     const {return m_cenzz; }

    /// Retrieve the centroid position
    inline const Vector getAxis()         const {return m_clusterAxis;}

    /// Errors in a HepMatrix
    const HepMatrix     getAxisErrs()     const;

    /// Direct access to errors
    inline const double getxDirxDir()     const {return m_axisxx; }
    inline const double getxDiryDir()     const {return m_axisxy; }
    inline const double getxDirzDir()     const {return m_axisxz; }
    inline const double getyDiryDir()     const {return m_axisyy; }
    inline const double getyDirzDir()     const {return m_axisyz; }
    inline const double getzDirzDir()     const {return m_axiszz; }

    /// Define an ( ) operator (allows read/write - indexing from 1!!)
    //double& operator()(const int &i);
    //const double& operator()(const int &i) const;
    //double& operator()(const int &i, const int &j);
    //const double& operator()(const int &i, const int &j) const;

    /// Set the energy
    inline void setEnergy(const double& energy)         {m_energy   = energy;}
    inline void setEnergyErr(const double& energyErr)   {m_eneError = energyErr;}

    /// Set parameters for centroid
    inline void setCentroid(const Point& pos)           {m_clusterCentroid = pos; }

    void        setCentroidErrs(const HepMatrix& errs);

    /// Set errors
    inline void setxDirxDir(const double& val)          {m_cenxx = val; }
    inline void setxDiryDir(const double& val)          {m_cenxy = val; }
    inline void setxDirzDir(const double& val)          {m_cenxz = val; }
    inline void setyDiryDir(const double& val)          {m_cenyy = val; }
    inline void setyDirzDir(const double& val)          {m_cenyz = val; }
    inline void setzDirzDir(const double& val)          {m_cenzz = val; }

    /// Set parameters for axis
    inline void setAxis(const Vector& axis)             {m_clusterAxis = axis; }

    void        setAxisErrs(const HepMatrix& errs);

    /// Set errors
    inline void setxPosxPos(const double& val)          {m_axisxx = val; }
    inline void setxPosyPos(const double& val)          {m_axisxy = val; }
    inline void setxPoszPos(const double& val)          {m_axisxz = val; }
    inline void setyPosySlp(const double& val)          {m_axisyy = val; }
    inline void setyPoszPos(const double& val)          {m_axisyz = val; }
    inline void setzPoszPos(const double& val)          {m_axiszz = val; }

    std::ostream& fillStream( std::ostream& s ) const;
    friend std::ostream& operator<< ( std::ostream& s, const CalParams& obj ) 
    {
	    return obj.fillStream(s);
    }

private:
    /// Private intialization method
    void initDataMembers();
  
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

