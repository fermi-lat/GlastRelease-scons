#ifndef test_util_h
#define test_util_h
// $Header$

// STD INCLUDES
#include <cmath>

/** @file 
    @author Z.Fewtrell
    generic utitilies for CalXtalResponse test apps.
*/

namespace CalXtalResponse {
  /// tolerance for reconstructed energy values
  /// \note - Due to asymmetry spline curvature in real calibrations, I have personally seen this effect as bad as 2%
  static const float MAX_RECON_DIFF = .01;

  /// since all ADC values are integers, they may be up to .5
  /// away from 'true' value.
  static const float MAX_ADC_DIFF = 1;

  /// slightly more forgiving for spline points
  static const float MAX_SPLINE_DIFF = .005;

  /// possible variation in asymmetry curves (muon based HE calibration has high error on each point)
  static const float MAX_ASYM_DIFF = .1;


  inline float abs_diff(const float a, const float b) {
    return std::fabs(a-b);
  }

  /// return relative diff (abs) between 2 floats avoid divide-by-zero
  /// errros unless
  inline float rel_diff(const float a, const float b) {
    // safe to divide by a
    if (a!=0) return std::fabs((a-b)/a);

    // fall back to divide by b
    if (b!=0) return std::fabs((a-b)/b);

    // only possibility a==b==0
    // zero pct diff
    return 0;
  }

  /// check that relative difference between a and b < tolerance
  /// if abs(a) < tolerance && abs(b) < tolerance, then return
  /// abs(a-b) < tolerance
  inline bool smart_compare(const float a,
                            const float b,
                            const float tolerance) {
    if (fabs(a) < tolerance && fabs(b) < tolerance)
      return fabs(a-b) < tolerance;

    else return rel_diff(a,b) < tolerance;
        
  }

  /// compare each value in 2 iterable collections for rel_diff() <
  /// tolerance, interface similar to STL <algorithm>
  template <typename _It>
  bool compare_collection(_It begin1,
                          _It end1,
                          _It begin2,
                          float tolerance) {
    _It it1(begin1);
    _It it2(begin2);
    while (it1 != end1) {
      if (rel_diff(*it1, *it2) > tolerance)
        return false;

      it1++, it2++;
    }

    return true;
  }

  /// length of csi crystal in mm
  static const float csiMM = 326;
        
  /// convert xtal (-167...167) position in mm to relative position (0...1)
  inline float xtalMM2RelPos(const float mm) {
    return mm/csiMM + .5;
  }

  /// convert xtal relative position (0...1) to mm (-167-167)
  inline float xtalRelPos2MM(const float relPos) {
    return (relPos -.5)*csiMM;
  }


  /// \brief ensure that a threshold test is correct (given a certain
  /// margin for error)
  /// \param test threshold level
  /// \param input signal
  /// \param measured test result
  /// \param margin threhold to ignore if signal is close to
  /// margin. (in same units as threshold)

  inline bool trig_test_margin(const float signal, const float thresh, const bool result, const double margin) {
    // return false if the trigger should have gone high, but it
    // didn't
    if (signal >= thresh + margin && !result) return false;
    // return true if it should have gone low but fired anyway
    if (signal < thresh - margin && result) return false;

    // otherwise return true
    return true;
  }
}

#endif // test_util_h
