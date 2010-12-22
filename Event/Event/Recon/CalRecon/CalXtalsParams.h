#ifndef CalXtalsParams_H
#define CalXtalsParams_H


/** 
* @class CalXtalsParams
*
* @brief Gaudi TDS class to store the basic quantities related to a collection of CAL xtals.
*
* These include the number of xtals in  the collection, the number of
* saturated xtals, the number ox xtals exceeding a certainf fraction of the total
* energy, the moments of the xtal energy distribution etc.
*
* Whenever this class is changed, the changes should be propagated to the
* related files on the ROOT side:
* - reconRootData/reconRootData/CalXtalsParams.h
* - reconRootData/src/CalXtalsParams.cxx
* 
* @author Luca Baldini, Johan Bregeon
*
*/

#include <iostream>


namespace Event { //Namespace Event
  
  
  class CalXtalsParams
  {
  public:
    /// Default (no parameter) constructor.
    CalXtalsParams() { clear(); }
    
    /// Constructor from all members.
    CalXtalsParams(int numXtals, int numTruncXtals, int numSaturatedXtals,
		   double xtalRawEneSum, double xtalCorrEneSum,
		   double xtalEneRms, double xtalEneSkewness);

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
    inline double getXtalEneRms()          const { return m_xtalEneRms; }
    inline double getXtalEneSkewness()     const { return m_xtalEneSkewness; }

    /// Set class parameters.
    inline void setNumXtals(int val)             { m_numXtals = val; }
    inline void setNumTruncXtals(int val)        { m_numTruncXtals = val; }
    inline void setNumSaturatedXtals(int val)    { m_numSaturatedXtals = val; }
    inline void setXtalRawEneSum(double val)     { m_xtalRawEneSum = val; }
    inline void setXtalCorrEneSum(double val)    { m_xtalCorrEneSum = val; }
    inline void setXtalEneRms(double val)        { m_xtalEneRms = val; }
    inline void setXtalEneSkewness(double val)   { m_xtalEneSkewness = val; }

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
    /// Rms of the xtal energy distribution.
    double m_xtalEneRms;
    /// Skewness of the xtal energy distribution.
    double m_xtalEneSkewness;
  };
  

}; //Namespace Event

#endif

