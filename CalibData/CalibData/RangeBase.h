// $Header$
#ifndef CalibData_RangeBase_h
#define CalibData_RangeBase_h
namespace CalibData {


  /** 
        Base class for per crystal-face-range Calorimeter calibration data
        or for per tile-pmt-range ACD calibration data
   */
  class RangeBase {

  public:
    RangeBase() {}
    virtual ~RangeBase() {}

    /// Derived classes will do a dynamic cast of argument, which 
    /// must be of same type, and then a deep copy.
    virtual void update(RangeBase* ) {}
    
  };

  /** 
       Generally speaking each value in a calorimeter calibration comes 
       with an associated uncertainty.  Put them together with this
       little class
  */

  class ValSig {
  public:
    ValSig(float val=-1, float sig=-1) : m_val(val), m_sig(sig) {} 
    ValSig(const ValSig& other) {m_val = other.m_val; m_sig = other.m_sig;}
    bool isDefined( ) const {return (m_sig >= 0.0); }
    float getVal() const {return m_val;}
    float getSig() const {return m_sig;}

    void setUndefined() {m_val= -1.0; m_sig = -1;}
    float m_val;
    float m_sig;
  };
}
#endif
