/**
 * @class CalMomentsAnalysis
 *
 * @brief Implements a "Moments Analysis" for use with categorizing tracker events before recon
 *        This is taken from code originally authored by Jay Norris and Heather Arrighi in 1998
 *        (see CalRecon for updated version of that code)and is based on the determination of an
 *        inertia tensor a la H. Goldstein in "Classical Mechanics", 1965. 
 *
 * @author Tracy Usher
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
    //@brief CalMomentsData is a utility data object for the moments analysis which 
    //       attempts to make the class independent of the actual Cal data objects used
    //       Minimum constructor requires position and weight for the data point
    CalMomentsData(const Point& point, const double weight, const double distToAxis=0.) :
                   m_useFlag(true), m_point(point), m_weight(weight), m_distToAxis(distToAxis) {};
    ~CalMomentsData() {}

    //@brief Provides access to data
    const Point&  getPoint()      const {return m_point;}
    const double  getWeight()     const {return m_weight;}
    const double  getDistToAxis() const {return m_distToAxis;}
    const double  getCoordAlongAxis() const {return m_coordAlongAxis;}
    bool          useIt()               {return m_useFlag;}

    //@brief Provides "set" functions
    void setPoint(const Point& point) {m_point   = point;}
    void setWeight(double weight)     {m_weight  = weight;}
    void setUseFlag(bool flag)        {m_useFlag = flag;}

    //@brief Determine distance to given axis
    double calcDistToAxis(const Point& centroid, const Vector& axis);

    //@brief Determine the coordinate along a given axis
    double calcCoordAlongAxis(const Point& centroid, const Vector& axis);

    // Define how to sort
    const bool operator<(const CalMomentsData& right)  const {return m_distToAxis < right.getDistToAxis();}

private:
    // bool for using or not using this data value
    bool   m_useFlag;
    // The position of this data point
    Point  m_point;
    // A weight to assign to the point in the moments calculation
    double m_weight;
    // The distance from the "axis" of this point
    double m_distToAxis;
    // The position along the "axis" of this point (with sign, used to calculate the skewness)
    double m_coordAlongAxis;
};

typedef std::vector<CalMomentsData> CalMomentsDataVec;

class CalMomentsAnalysis
{
public:
    //@brief Define the CalMomentsAnalysis class
    /// Constructor takes no arguments
    CalMomentsAnalysis();
    ~CalMomentsAnalysis() {};

    //@brief Perform the moments analysis on the given data around the given point
    double doMomentsAnalysis(CalMomentsDataVec& dataVec, const Point& centroid);

    //@brief Drives an iterative moments analysis
    //       Note the input data vector is NOT a reference (so is a copy)
    double doIterativeMomentsAnalysis(CalMomentsDataVec dataVec, 
                                      const Point&      centroid,
                                      double            scaleFactor);

    //@brief Extract the results
    const Point  getMomentsCentroid()       const {return m_centroid;}
    const Vector getMoments()               const {return m_moment;}
    const Vector getMomentsAxis(int axis=1) const {return m_axis[axis];}
    const double getLongitudinalRms()       const {return m_rmsLong;}
    const double getTransverseRms()         const {return m_rmsTrans;}
    const double getLongAsymmetry()         const {return m_rmsLongAsym;}
    const double getWeightSum()             const {return m_weightSum;}
    const double getLongSkewness()          const {return m_skewnessLong;}
    const double getNumIterations()         const {return m_numIterations;}
    const double getNumDroppedPoints()      const {return m_numDroppedPoints;}

private:

    // Centroid of the moments
    Point  m_centroid;
    // Vector of calculated moments
    Vector m_moment;
    // Axis corresponding to the longest principal moment
    Vector m_axis[3];

    // The Longitudinal rms of the distribution of points
    double m_rmsLong;
    // The transverse rms of the distribution of points
    double m_rmsTrans;
    // The asymmetry associated with the longitudinal distribution
    double m_rmsLongAsym;
    // Sum of weights in moments analysis 
    double m_weightSum;
    // Skewness in the moments analysis
    double m_skewnessLong;
    // Statistics on iterations (if done)
    int    m_numIterations;
    int    m_numDroppedPoints;
};

