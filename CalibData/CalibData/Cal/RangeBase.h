// $Header$

namespace CalibData {

  /** Base class for per crystal-face-range Calorimeter calibration data.
   */
  class RangeBase {

  public:
    RangeBase() {}
    virtual ~RangeBase() {}

    /// Derived classes will do a dynamic cast of argument, which 
    /// must be of same type, and then a deep copy.
    virtual void update(RangeBase* ) {}
    
  };

}
