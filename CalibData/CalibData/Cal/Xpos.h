// $Header$
#ifndef CalibData_Xpos_h
#define CalibData_Xpos_h
#include <vector>
namespace CalibData {

  /**   @class Xpos
      Keeps track of a collection of x positions along a crystal. Normally
      a single Xpos will apply to all crystals being calibrated.
   */
  class Xpos {

  public:
    Xpos(std::vector<float>* vals=0);
    ~Xpos() {}

    const std::vector<float>* getVals() const {return &m_vals;}
    /// Deep copy.
    void update(const Xpos* other);

  private:
    std::vector<float> m_vals;
    
  };

}
#endif
