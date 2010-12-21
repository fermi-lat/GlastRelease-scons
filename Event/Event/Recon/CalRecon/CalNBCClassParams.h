#ifndef CalNBCClassParams_H
#define CalNBCClassParams_H


/** 
* @class CalNBCClassParams
*
* @brief Gaudi TDS class to store the output of the Naive Bayes Classifier.
*
* The class encapsulates the calculation and normalization of the probabilities for
* the NBC.
* 
* Whenever this class is changed, the changes should be propagated to the
* related files on the ROOT side:
* - reconRootData/reconRootData/CalNBCClassParams.h
* - reconRootData/src/CalNBCClassParams.cxx
* 
* @author Luca Baldini, Johan Bregeon
*
*/

#include <iostream>
#include "CalClassParams.h"


namespace Event { //Namespace Event
  
  
  class CalNBCClassParams : public CalClassParams
  {
  public:
    /// Default (no parameter) constructor.
    CalNBCClassParams();

    /// Constructor from all members.
    CalNBCClassParams(std::string producerName, std::map <std::string, double> probMap);

    /// Destructor.
    ~CalNBCClassParams() {}

    /// Handle the NBC algorithm basic steps.
    /// Multiply the probability value of the given class name by the pdf value from the xml
    /// (or set it to the pdf value dealing with the first variable in the list).
    void multiply(const std::string &className, double pdfValue);
    /// Normalize to 1 the sum of the probability values.
    void normalize();

  private:


  };


}; //Namespace Event

#endif
