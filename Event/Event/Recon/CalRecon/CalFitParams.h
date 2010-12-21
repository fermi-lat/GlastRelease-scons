
#ifndef CalFitParams_H
#define CalFitParams_H

#include "geometry/Point.h"
#include "geometry/Vector.h"
#include "CLHEP/Matrix/Matrix.h"

/** 
* @class CalFitParams
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


class CalFitParams
{
public:
    /// Default constructor
    CalFitParams() { clear() ; }

    /// reset method
    void clear() ;
  
    /// Constructor given Points, Vectors and matrices
    CalFitParams(int           nFitLayers, double                  chiSquare,
                 const Point&  centroid,   const CLHEP::HepMatrix& centroidErr,
                 const Vector& axis,       const CLHEP::HepMatrix& axisErr);

    /// Direct construction from all the elements (the old fashioned way)
    CalFitParams(int    nFitLayers, double chiSquare,
                 double xCntrd, double yCntrd, double zCntrd,
                 double cntdxx, double cntdxy, double cntdxz, double cntdyy, double cntdyz, double cntdzz,
                 double xAxis,  double yAxis,  double zAxis,
                 double axsdxx, double axsdxy, double axsdxz, double axsdyy, double axsdyz, double axsdzz);

    /// Copy constructor
    //CalFitParams (const CalFitParams& right);

   ~CalFitParams() {}

    /// Retrieve the energy
    inline int           getFitLayers()    const {return m_nFitLayers;}
    inline double        getChiSquare()    const {return m_chiSquare;}

    /// Retrieve the centroid position
    inline const Point&  getCentroid()     const {return m_clusterCentroid;}

    /// Errors in a HepMatrix
    CLHEP::HepMatrix     getCentroidErrs() const ;

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
    CLHEP::HepMatrix     getAxisErrs()     const ;

    /// Direct access to errors
    inline double        getxDirxDir()     const {return m_axisxx; }
    inline double        getxDiryDir()     const {return m_axisxy; }
    inline double        getxDirzDir()     const {return m_axisxz; }
    inline double        getyDiryDir()     const {return m_axisyy; }
    inline double        getyDirzDir()     const {return m_axisyz; }
    inline double        getzDirzDir()     const {return m_axiszz; }

    /// Set the parameters
    inline void setFitLayers(int nFitLayers)    {m_nFitLayers = nFitLayers;}
    inline void setChiSquare(double chiSquare)  {m_chiSquare  = chiSquare;}

    /// Set parameters for centroid
    inline void setCentroid(const Point& pos)   {m_clusterCentroid = pos; }

    void        setCentroidErrs(const CLHEP::HepMatrix& errs);

    /// Set errors
    inline void setxDirxDir(double val)         {m_cenxx = val; }
    inline void setxDiryDir(double val)         {m_cenxy = val; }
    inline void setxDirzDir(double val)         {m_cenxz = val; }
    inline void setyDiryDir(double val)         {m_cenyy = val; }
    inline void setyDirzDir(double val)         {m_cenyz = val; }
    inline void setzDirzDir(double val)         {m_cenzz = val; }

    /// Set parameters for axis
    inline void setAxis(const Vector& axis)     {m_clusterAxis = axis; }

    void        setAxisErrs(const CLHEP::HepMatrix& errs);

    /// Set errors
    inline void setxPosxPos(double val)         {m_axisxx = val; }
    inline void setxPosyPos(double val)         {m_axisxy = val; }
    inline void setxPoszPos(double val)         {m_axisxz = val; }
    inline void setyPosySlp(double val)         {m_axisyy = val; }
    inline void setyPoszPos(double val)         {m_axisyz = val; }
    inline void setzPoszPos(double val)         {m_axiszz = val; }

    std::ostream& fillStream( std::ostream& s ) const;
    friend std::ostream& operator<< ( std::ostream& s, const CalFitParams& obj ) 
    {
	    return obj.fillStream(s);
    }

private:

    /// Parameters associated with the fit
    int    m_nFitLayers;
    double m_chiSquare;

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


inline CalFitParams::CalFitParams(int   fitLayers, double chiSquare,
                                  const Point&  centroid, const CLHEP::HepMatrix& centroidErr,
                                  const Vector& axis,     const CLHEP::HepMatrix& axisErr) :
                            m_nFitLayers(fitLayers),
                            m_chiSquare(chiSquare),
                            m_clusterCentroid(centroid),
                            m_clusterAxis(axis)
{
    setCentroidErrs(centroidErr);
    setAxisErrs(axisErr);
}

inline CalFitParams::CalFitParams(int    fitLayers, double chiSquare,
                                  double xCntrd, double yCntrd, double zCntrd,
                                  double cntdxx, double cntdxy, double cntdxz, double cntdyy, double cntdyz, double cntdzz,
                                  double xAxis,  double yAxis,  double zAxis,
                                  double axsdxx, double axsdxy, double axsdxz, double axsdyy, double axsdyz, double axsdzz) :
                                  m_nFitLayers(fitLayers),
                                  m_chiSquare(chiSquare),
                                  m_cenxx(cntdxx),
                                  m_cenxy(cntdxy),
                                  m_cenxz(cntdxz),
                                  m_cenyy(cntdyy),
                                  m_cenyz(cntdyz),
                                  m_cenzz(cntdzz),
                                  m_axisxx(axsdxx),
                                  m_axisxy(axsdxy),
                                  m_axisxz(axsdxz),
                                  m_axisyy(axsdyy),
                                  m_axisyz(axsdyz),
                                  m_axiszz(axsdzz)
{
    m_clusterCentroid = Point(xCntrd,yCntrd,zCntrd);
    m_clusterAxis     = Vector(xAxis,yAxis,zAxis);
    return;
}

inline void CalFitParams::clear()
{
    m_clusterCentroid = Point(0.,0.,0.);
    m_clusterAxis     = Vector(0.,0.,0.);

    m_cenxx      = 0.;
    m_cenxy      = 0.;
    m_cenxz      = 0.;
    m_cenyy      = 0.;
    m_cenyz      = 0.;
    m_cenzz      = 0.;

    m_axisxx     = 0.;
    m_axisxy     = 0.;
    m_axisxz     = 0.;
    m_axisyy     = 0.;
    m_axisyz     = 0.;
    m_axiszz     = 0.;

    m_nFitLayers = 0;
    m_chiSquare  = 0.;
}

inline void CalFitParams::setCentroidErrs(const CLHEP::HepMatrix& errs)
{
    m_cenxx = errs(1,1);
    m_cenxy = errs(1,2);
    m_cenxz = errs(1,3);
    m_cenyy = errs(2,2);
    m_cenyz = errs(2,3);
    m_cenzz = errs(3,3);
}

inline void CalFitParams::setAxisErrs(const CLHEP::HepMatrix& errs)
{
    m_axisxx = errs(1,1);
    m_axisxy = errs(1,2);
    m_axisxz = errs(1,3);
    m_axisyy = errs(2,2);
    m_axisyz = errs(2,3);
    m_axiszz = errs(3,3);
}

inline CLHEP::HepMatrix CalFitParams::getCentroidErrs() const
{
    CLHEP::HepMatrix errs(3,3,0);

    errs(1,1) = m_cenxx;
    errs(1,2) = errs(2,1) = m_cenxy;
    errs(1,3) = errs(3,1) = m_cenxz;
    errs(2,2) = m_cenyy;
    errs(2,3) = errs(3,2) = m_cenyz;
    errs(3,3) = m_cenzz;

    return errs;
}

inline CLHEP::HepMatrix CalFitParams::getAxisErrs() const
{
    CLHEP::HepMatrix errs(3,3,0);

    errs(1,1) = m_axisxx;
    errs(1,2) = errs(2,1) = m_axisxy;
    errs(1,3) = errs(3,1) = m_axisxz;
    errs(2,2) = m_axisyy;
    errs(2,3) = errs(3,2) = m_axisyz;
    errs(3,3) = m_axiszz;

    return errs;
}

inline std::ostream& CalFitParams::fillStream( std::ostream& s ) const
{
  s <<
    "Centroid = (" << m_clusterCentroid.x() << ", " << m_clusterCentroid.y() << ", "
		<< m_clusterCentroid.z() << ") mm\n" <<
    "Centroid covariance matrix:\n" <<
    "| " << m_cenxx  << "  " << m_cenxy << "  " << m_cenxz << " |\n" <<
    "| " << m_cenxy  << "  " << m_cenyy << "  " << m_cenyz << " |\n" <<
    "| " << m_cenxz  << "  " << m_cenyz << "  " << m_cenzz << " |\n" <<
    "Axis = (" << m_clusterAxis.x() << ", " << m_clusterAxis.y() << ", "
		  << m_clusterAxis.z() << ")\n" <<
    "Axis covariance matrix:\n" <<
    "| " << m_axisxx  << "  " << m_axisxy << "  " << m_axisxz << " |\n" <<
    "| " << m_axisxy  << "  " << m_axisyy << "  " << m_axisyz << " |\n" <<
    "| " << m_axisxz  << "  " << m_axisyz << "  " << m_axiszz << " |\n" <<
    "Number of layers for the fit: " << m_nFitLayers << "\n" <<
    "Fit chisquare: " << m_chiSquare;

    return s;
}


}; //Namespace Event


#endif

