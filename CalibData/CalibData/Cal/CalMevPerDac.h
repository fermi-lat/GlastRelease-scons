// $Header$

/// @file CalMevPerDac.h
/** 
  This file contains class definitions corresponding to the calibration
  type  CAL_MevPerDac.    Two classes are defined.  One, CalMevPerDac, 
  is the data pertaining to a  single crystal. The other, CalMevPerDacCol, 
  is a collection of CalMevPerDac instances. 

 @author J. Bogart
*/

#ifndef CalibData_CalMevPerDac_h
#define CalibData_CalMevPerDac_h

#include "CalibData/RangeBase.h"
#include "CalibData/Cal/CalCalibBase.h"

namespace CalibData {
  /** @class CalMevPerDac
     Data pertaining to a single crystal.  It inherits from RangeBase.
  */
  class CalMevPerDac : public RangeBase {
  public:
    CalMevPerDac(const ValSig* big = 0, 
               const ValSig* small = 0,
               const std::vector<ValSig>* bigSmallRatioN = 0,
               const std::vector<ValSig>* bigSmallRatioP = 0);
    ~CalMevPerDac();

    // Bunch of get functions

    const ValSig* getBig() const {return &m_big;}
    const ValSig* getSmall() const {return &m_small;}

    /// Caller should use enum definitions in idents::CalXtalId::XtalFace
    /// for @a face argument.
    const std::vector<ValSig>* getBigSmallRatio(int face) const;

    virtual void update(RangeBase* other);

  private:

    /// Gain and uncertainty for sqrt (P*N), big diode
    ValSig m_big;
    /// Gain and uncertainty for sqrt (P*N), small diode
    ValSig m_small;

    /// Ratio of big diode/small diode for Neg end as a function of position
    std::vector<ValSig>* m_bigSmallRatioN;
    /// Ratio of big diode/small diode for Pos end as a function of position
    std::vector<ValSig>* m_bigSmallRatioP;

  };

  /**  @class CalMevPerDacCol
     This class represents CalMevPerDac information for the entire detector.
     Individual CalMevPerDac instances may be looked up by CalXtalId.
  */
  class CalMevPerDacCol :  public CalCalibBase {

  public:
    CalMevPerDacCol(unsigned nTowerRow=4, unsigned nTowerCol=4, 
                    unsigned nLayer=8, unsigned nXtal=12, unsigned nXpos=1);

    ~CalMevPerDacCol();

    /**
       Store data for a single channel.  Most of the work is done
       by the base class @a CalCalibBase
    */
    bool putRange(idents::CalXtalId id, unsigned range, 
                  unsigned face, RangeBase* data);

    /**
       Store data for a single channel.  Most of the work is done
       by the base class @a CalCalibBase
    */
    bool putRange(unsigned towerRow, unsigned towerCol, 
                  unsigned layer, unsigned xtal, unsigned range,
                  unsigned face, RangeBase* data);

    virtual const CLID& clID() const {return classID(); }

    static const CLID& classID();

  };
}
#endif
