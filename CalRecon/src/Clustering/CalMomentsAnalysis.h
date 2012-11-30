/**
 * @class CalMomentsAnalysis
 *
 * @brief Implements a "Moments Analysis" for use with categorizing tracker events before recon
 *        This is taken from code originally authored by Jay Norris and Heather Arrighi in 1998
 *        (see CalRecon for updated version of that code)and is based on the determination of an
 *        inertia tensor a la H. Goldstein in "Classical Mechanics", 1965. 
 *
 * @author Tracy Usher, Luca Baldini.
 *
 *
 * File and Version Information:
 *      $Header$
 */

#include "geometry/Ray.h"
#include "Event/Recon/CalRecon/CalXtalRecData.h"
#include "CalMomentsData.h"

#include <vector>

class CalMomentsAnalysis
{
 public:
  /// No-parameter constructor
  CalMomentsAnalysis() { clear(); }
  
  /// Destructor
  ~CalMomentsAnalysis() {};

  /// Reset.
  void clear();
  
  /// Perform the moments analysis on the given data around the given point
  double doMomentsAnalysis(CalMomentsDataVec& dataVec, const Point& centroid,
                           double coreRadius);

  /// Drives an iterative moments analysis
  /// (note the input data vector is NOT a reference, it is a copy)
  double doIterativeMomentsAnalysis(CalMomentsDataVec dataVec, const Point& centroid,
                                    double transScaleFactor, double transScaleFactorBoost,
                                    double coreRadius);

  /// ADW: Sept. 7, 2011
  /// Calculate the covariance on the axis direction
  /// CS: Nov 2012, update to an implementation that uses xtal errors
  void calcCovarianceAxisSimple(Vector momAxis);
  void calcCovarianceCentroidSimple(CalMomentsDataVec& dataVec, Vector momAxis, const Point& centroid);
  void calcCovariance(CalMomentsDataVec& dataVec, const Point& centroid);

  /// Access class members...
  inline const double getWeightSum()        const { return m_weightSum; }
  inline const Point  getCentroid()         const { return m_centroid; }
  inline const CLHEP::HepMatrix  getCentroidErr()  const { return m_centroidErr; }
  inline const Vector getMoments()          const { return m_moment; }
  inline const Vector getAxis(int axis=1)   const { return m_axis[axis]; }
  inline const CLHEP::HepMatrix  getAxisErr()  const { return m_axisErr; }
  inline const double getFullLength()       const { return m_fullLength; }
  inline const double getLongRms()          const { return m_longRms; }
  inline const double getTransRms()         const { return m_transRms; }
  inline const double getLongRmsAsym()      const { return m_longRmsAsym; }
  inline const double getLongSkewness()     const { return m_longSkewness; }
  inline const double getCoreEnergyFrac()   const { return m_coreEnergyFrac; }
  inline const int    getNumIterations()    const { return m_numIterations; }
  inline const int    getNumDroppedPoints() const { return m_numDroppedPoints; }

 private:
  
  /// Sum of weights in moments analysis 
  double m_weightSum;
  /// Centroid of the cluster.
  Point m_centroid;
  /// ADW: Sept. 7, 2011
  /// Matrix corresponding to the error on the centoid
  CLHEP::HepMatrix m_centroidErr;
  /// Vector of calculated moments.
  Vector m_moment;
  /// Axis corresponding to the principal moments.
  Vector m_axis[3];
  /// ADW: Sept. 7, 2011
  /// Matrix correspoding to the covariance of the principal axis
  CLHEP::HepMatrix m_axisErr;
  /// The distance (in radiation lengths) between the positions of the first and
  /// the last xtal, projected along the main axis of the cluster.
  double m_fullLength;
  /// The Longitudinal rms of the distribution of points.
  double m_longRms;
  /// The transverse rms of the distribution of points.
  double m_transRms;
  /// The asymmetry associated with the longitudinal distribution
  double m_longRmsAsym;
  /// Skewness in the moments analysis
  double m_longSkewness;
  /// Fractional energy sum within 1 (or wathever) Moliere radius from the cluster axis.
  double m_coreEnergyFrac;
  /// Statistics on iterations (if done)
  int m_numIterations;
  /// Number of calorimeter hits dropped during the iterations.
  int m_numDroppedPoints;
};

