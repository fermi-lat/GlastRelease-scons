
#ifndef CalMomParams_H
#define CalMomParams_H


/** 
* @class CalMomParams
*
* @brief Gaudi TDS class to store the complete output of the calorimeter
* moments analysis.
*
* This is an extension of the CalParams base class meant to include all the
* quantities calculated in the moments analysis (longitudinal and transverse
* RMS of the cluster, skewness, moments asymmetry and so on and so forth)
* along with the basic information on the cluster energy, centroid and
* direction.
*
* Whenever this class is changed, the changes should be propagated to the
* related files on the ROOT side:
* - reconRootData/reconRootData/CalMomParams.h
* - reconRootData/src/CalMomParams.cxx
* 
* @author Luca Baldini, Johan Bregeon
*
*/

#include <iostream>
#include "CalParams.h"


namespace Event { //Namespace Event
  

  class CalMomParams : public CalParams
  {
  public:
    /// Default (no parameter) constructor.
    CalMomParams() { clear() ; }

    /// Direct construction from all the elements.
    CalMomParams(double energy, double eneError,
		 const Point& centroid, const CLHEP::HepMatrix& centroidErr,
		 const Vector& axis, const CLHEP::HepMatrix& axisErr,
		 int numIterations, double transRms, double longRms, double longRmsAsym,
		 double longSkewness, double coreEnergyFrac);

    /// And even more parameters (reflecting the old-fashioned way CalParams constructor).
    CalMomParams(double energy, double eneError,
		 double xCntrd, double yCntrd, double zCntrd,
		 double cntdxx, double cntdxy, double cntdxz,
		 double cntdyy, double cntdyz, double cntdzz,
		 double xAxis,  double yAxis,  double zAxis,
		 double axsdxx, double axsdxy, double axsdxz,
		 double axsdyy, double axsdyz, double axsdzz,
		 int numIterations, double transRms, double longRms, double longRmsAsym,
		 double longSkewness, double coreEnergyFrac);

    /// Convenience constructor to be used to replace an old CalParams object directly
    /// (i.e. the specific CalMomParams members are automagically initialized).
    CalMomParams(double energy, double eneError,
		 double xCntrd, double yCntrd, double zCntrd,
		 double cntdxx, double cntdxy, double cntdxz,
		 double cntdyy, double cntdyz, double cntdzz,
		 double xAxis,  double yAxis,  double zAxis,
		 double axsdxx, double axsdxy, double axsdxz,
		 double axsdyy, double axsdyz, double axsdzz);

    /// Destructor.
    ~CalMomParams() {}
    
    /// Reset method.
    void clear();
    /// Part of the reset code specific to the CalMomParams class
    /// (as opposed to the base class CalParams).
    void clearMomParams();

    /// Retrieve parameters...
    inline double getTransRms()       const {return m_transRms;}
    inline double getLongRms()        const {return m_longRms;}
    inline double getLongRmsAsym()    const {return m_longRmsAsym;}
    inline double getLongSkewness()   const {return m_longSkewness;}
    inline double getCoreEnergyFrac() const {return m_coreEnergyFrac;}

    /// Return the ratio between the tranverse and the longitudinal RMS values.
    double getElongation() const;

    /// Set parameters.
    inline void setTransRms(double val)       {m_transRms = val;}
    inline void setLongRms(double val)        {m_longRms = val;}
    inline void setLongRmsAsym(double val)    {m_longRmsAsym = val;}
    inline void setLongSkewness(double val)   {m_longSkewness = val;}
    inline void setCoreEnergyFrac(double val) {m_coreEnergyFrac = val;}

    /// Std output facility.
    std::ostream& fillStream(std::ostream& s) const;
    friend std::ostream& operator<< (std::ostream& s, const CalMomParams& obj)
    {
      return obj.fillStream(s);
    }

  private:
    
    /// The number of iterations in the moment analysis.
    int m_numIterations;
    /// TBD add more statistics on the iterations (i.e. the number of xtals at the last step).
    /// The transverse RMS of the energy distribution in the cluster.
    double m_transRms;
    /// The longitudinal RMS of the energy distribution in the cluster
    /// (the average of the two largest moments of the distribution).
    double m_longRms;
    /// The longitudinal RMS asymmetry (the fractional difference between
    /// the two largest moments of the distribution).
    double m_longRmsAsym;
    /// The skewness of the energy profile along the cluster axis.
    double m_longSkewness;
    /// Fractional energy sum within 1 Moliere radius from the cluster axis.
    double m_coreEnergyFrac;
  };


}; //Namespace Event

#endif
