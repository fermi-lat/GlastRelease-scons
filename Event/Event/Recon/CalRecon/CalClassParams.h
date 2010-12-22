#ifndef CalClassParams_H
#define CalClassParams_H


/** 
* @class CalClassParams
*
* @brief Gaudi TDS class to store the output of the cluster classification.
*
* The basic member is a std::map <std:string, double> containing the probability
* value for each topological class, as return by any generic classification algorithm.
* 
* The list of classes is a priori undetermined (it's up to the classifier),
* except for the fact that there's *always* a "gam" class whose probability is
* initialized to -1. in the constructor.
*
* Whenever this class is changed, the changes should be propagated to the
* related files on the ROOT side:
* - reconRootData/reconRootData/CalClassParams.h
* - reconRootData/src/CalClassParams.cxx
* 
* @author Luca Baldini, Johan Bregeon
*
*/

#include <iostream>
#include <map>


namespace Event { //Namespace Event
  
  
  class CalClassParams
  {
  public:
    /// Default (no parameter) constructor.
    CalClassParams();

    /// Constructor from all members.
    CalClassParams(std::string producerName, std::map <std::string, double> probMap);

    /// Destructor.
    ~CalClassParams() {}
    
    /// Reset method.
    void clear();

    /// Retrieve class parameters...
    inline std::map <std::string, double> getProbMap() const { return m_probMap; }
    inline std::string getProducerName()               const { return m_producerName; }
    double getProb(const std::string &className)       const;
    double getGamProb()                                const { return getProb("gam"); }

    /// Set class parameters.
    inline void setProbMap(std::map <std::string, double> probMap) { m_probMap = probMap; }
    inline void setProducerName(std::string producerName)    { m_producerName = producerName; }
    inline void setProb(std::string className, double prob)  { m_probMap[className] = prob; }

    /// Check whether the base std::map container has already a given key.
    bool hasClass(const std::string &className) const;

    /// Std output facility.
    std::ostream& fillStream(std::ostream& s) const;
    friend std::ostream& operator<< (std::ostream& s, const CalClassParams& obj)
    {
      return obj.fillStream(s);
    }

  private:

    /// Name of the producer.
    std::string m_producerName;
    /// The std::map containing the probability values for the different topologies.
    std::map <std::string, double> m_probMap;
  };


}; //Namespace Event

#endif
