#ifndef CalXtalsParams_H
#define CalXtalsParams_H


/**
   @file CalXtalsParams.h
   @class CalXtalsParams

   @brief Gaudi TDS class to store the basic quantities related to a collection of CAL xtals.

   These include the number of xtals in  the collection, the number of
   saturated xtals, the number ox xtals exceeding a certainf fraction of the total
   energy, the moments of the xtal energy distribution etc.
   
   Whenever this class is changed, the changes should be propagated to the
   related files on the ROOT side:
   - reconRootData/reconRootData/CalXtalsParams.h
   - reconRootData/src/CalXtalsParams.cxx
   
   @author Luca Baldini (luca.baldini@pi.infn.it)

   $Revision$
   $Date$
   $Header$
*/

#include <iostream>
#include "geometry/Point.h"


namespace Event { //Namespace Event
  
  
  class CalXtalsParams
  {
  public:
    /// Default (no parameter) constructor.
    CalXtalsParams() { clear(); }
    
    /// Constructor from all members.
    CalXtalsParams(int numXtals, int numTruncXtals, int numSaturatedXtals,
		   double xtalRawEneSum, double xtalCorrEneSum, double xtalEneMax,
		   double xtalEneRms, double xtalEneSkewness, const Point& centroid);

    /// Convenience constructor from the two members originally in the CalCluster class.
    CalXtalsParams(int numTruncXtals, int numSaturatedXtals);
    
    /// Destructor.
    ~CalXtalsParams() {}
    
    /// Reset method.
    void clear();

    /// Retrieve class parameters...
    inline int getNumXtals()               const { return m_numXtals; }
    inline int getNumTruncXtals()          const { return m_numTruncXtals; }
    inline int getNumSaturatedXtals()      const { return m_numSaturatedXtals; }
    inline double getXtalRawEneSum()       const { return m_xtalRawEneSum; }
    inline double getXtalCorrEneSum()      const { return m_xtalCorrEneSum; }
    inline double getXtalEneMax()          const { return m_xtalEneMax; }
    inline double getXtalEneRms()          const { return m_xtalEneRms; }
    inline double getXtalEneSkewness()     const { return m_xtalEneSkewness; }
    inline const Point&  getCentroid()     const { return m_centroid;}

    /// Set class parameters.
    inline void setNumXtals(int val)             { m_numXtals = val; }
    inline void setNumTruncXtals(int val)        { m_numTruncXtals = val; }
    inline void setNumSaturatedXtals(int val)    { m_numSaturatedXtals = val; }
    inline void setXtalRawEneSum(double val)     { m_xtalRawEneSum = val; }
    inline void setXtalCorrEneSum(double val)    { m_xtalCorrEneSum = val; }
    inline void setXtalEneMax(double val)        { m_xtalEneMax = val; }
    inline void setXtalEneRms(double val)        { m_xtalEneRms = val; }
    inline void setXtalEneSkewness(double val)   { m_xtalEneSkewness = val; }
    inline void setCentroid(const Point& pos)    { m_centroid = pos; }

    /// Std output facility.
    std::ostream& fillStream(std::ostream& s) const;
    friend std::ostream& operator<< (std::ostream& s, const CalXtalsParams& obj)
    {
      return obj.fillStream(s);
    }

  private:

    /// Number of xtals.
    int m_numXtals;
    /// Number of Xtals with > 1% (adjustable) of the total cluster energy.
    int m_numTruncXtals;
    /// Number of saturated xtals.
    int m_numSaturatedXtals;
    /// Plain sum of the xtal energies.
    double m_xtalRawEneSum;
    /// Corrected sum of the xtal energies.
    double m_xtalCorrEneSum;
    /// The energy of the xtal with the highest signal.
    double m_xtalEneMax;
    /// Rms of the xtal energy distribution.
    double m_xtalEneRms;
    /// Skewness of the xtal energy distribution.
    double m_xtalEneSkewness;
    /// Centroid of the xtal collection.
    Point m_centroid;
  };
  

}; //Namespace Event

#endif

