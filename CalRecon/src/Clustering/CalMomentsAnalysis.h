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
 * File and Version Information:
 *      $Header$
 */

#include "geometry/Ray.h"

#include <vector>

// Start by defining a data class to operate on
class CalMomentsData
{
 public:
  /// CalMomentsData is a utility data object for the moments analysis which 
  /// attempts to make the class independent of the actual Cal data objects used.
  /// Minimum constructor requires position and weight for the data point.
  CalMomentsData(const Point& point, const double weight, int tower) :
  m_useFlag(true),
    m_point(point),
    m_weight(weight),
    m_tower(tower),
    m_distToAxis(0.),
    m_coordAlongAxis(0.)
    {};
    
  /// Destructor
  ~CalMomentsData() {}
    
  /// Provides access to data.
  bool useIt()                     const { return m_useFlag; }
  const Point& getPoint()          const { return m_point; }
  const double getWeight()         const { return m_weight; }
  const int getTower()             const { return m_tower; }
  const double getDistToAxis()     const { return m_distToAxis; }
  const double getCoordAlongAxis() const { return m_coordAlongAxis; }
  
  /// Provides "set" functions.
  void setPoint(const Point& point)      { m_point = point; }
  void setWeight(double weight)          { m_weight = weight; }
  void setUseFlag(bool flag)             { m_useFlag = flag; }
  void setTower(int tower)               { m_tower = tower; }
  void setDistToAxis(double dist)        { m_distToAxis = dist; }
  void setCoordAlongAxis(double coord)   { m_coordAlongAxis = coord; }
  
  /// Determine distance to given axis.
  /// Note that this sets a class member, as well.
  double calcDistToAxis(const Point& centroid, const Vector& axis);
  
  /// Determine the coordinate along a given axis.
  /// Note that this sets a class member, as well.
  double calcCoordAlongAxis(const Point& centroid, const Vector& axis);
  
  /// Define how to sort.
  const bool operator<(const CalMomentsData& right)  const
  {return m_distToAxis < right.getDistToAxis();}
  
 private:
  /// Bool for using or not using this data value.
  bool m_useFlag;
  /// The position of this data point.
  Point m_point;
  /// A weight to assign to the point in the moments calculation.
  double m_weight;
  /// The tower this hit belongs to. This van be used in the CAL moments analysis
  /// to figure out whether two hits belong to the same tower or not.
  int m_tower;
  /// The distance from the "axis" of this point.
  double m_distToAxis;
  /// The position along the "axis" of this point (with sign, used to calculate the skewness).
  double m_coordAlongAxis;
};

typedef std::vector<CalMomentsData> CalMomentsDataVec;


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
  
  /// Access class members...
  inline const double getWeightSum()        const { return m_weightSum; }
  inline const Point  getCentroid()         const { return m_centroid; }
  inline const Vector getMoments()          const { return m_moment; }
  inline const Vector getAxis(int axis=1)   const { return m_axis[axis]; }
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
  /// Vector of calculated moments.
  Vector m_moment;
  /// Axis corresponding to the principal moments.
  Vector m_axis[3];
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

