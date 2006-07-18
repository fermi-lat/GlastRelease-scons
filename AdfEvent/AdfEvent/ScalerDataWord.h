#ifndef SCALERDATAWORD_HH
#define SCALERDATAWORD_HH

#include <iostream>
#include "AdfEvent/TaggerParameters.h"
/**
* @class ScalerDataWord 
* @brief Base class describing a data word from the scaler.

* The basic structure is as follows:
* 31          0
* -------------
* counter_value
*/

namespace AncillaryData {
  class ScalerDataWord {
  public:
    ScalerDataWord()                 {;} 
    ~ScalerDataWord()                 {;} 
    void setData(unsigned int word) {m_word = word;}
    unsigned int getScalerValue()    {return m_word;}
  private:
    unsigned int m_word;
  };
  
}//namespace AdfEvent
#endif
