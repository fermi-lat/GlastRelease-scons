// $Header$

/// @file CalAsym.h
/** 
  This file contains class definitions corresponding to the calibration
  type  CAL_Asym.    Two classes are defined.  One, CalAsym, 
  is the data pertaining to a  single crystal. The other, CalAsymCol, 
  is a collection of CalAsym instances. 

 @author J. Bogart
*/

#ifndef CalibData_CalAsym_h
#define CalibData_CalAsym_h

#include "CalibData/RangeBase.h"
#include "CalibData/Cal/CalCalibBase.h"

namespace CalibData {
  /** @class CalAsym
     Data pertaining to a single crystal.  It inherits from RangeBase.
  */
  class CalAsym : public RangeBase {
  public:
    CalAsym(const std::vector<ValSig>* big = 0,
            const std::vector<ValSig>* small = 0,
            const std::vector<ValSig>* nSmallPBig = 0,
            const std::vector<ValSig>* pSmallNBig = 0);
    ~CalAsym();

    const std::vector<ValSig>* getBig() const {return m_big;}
    const std::vector<ValSig>* getSmall() const {return m_small;}
    const std::vector<ValSig>* getNSmallPBig() const {return m_nSmallPBig;}
    const std::vector<ValSig>* getPSmallNBig() const {return m_pSmallNBig;}

    virtual void update(RangeBase* other);

  private:
    /// signal asymmetry over position for big diodes on both faces
    std::vector<ValSig>* m_big;
    /// signal asymmetry over position for small diodes on both faces
    std::vector<ValSig>* m_small;
    /// signal asymm over position w/ pos face big and minus face small diodes
    std::vector<ValSig>* m_nSmallPBig;
    /// signal asymm over position w/ pos face small and minus face big diodes
    std::vector<ValSig>* m_pSmallNBig;
  };

  /**  @class CalAsymCol
     This class represents CalAsym information for the entire detector.
     Individual CalAsym instances may be looked up by CalXtalId.
  */
  class CalAsymCol :  public CalCalibBase {

  public:
    CalAsymCol(unsigned nTowerRow=4, unsigned nTowerCol=4, 
               unsigned nLayer=8, unsigned nXtal=12, unsigned nXpos=1);

    ~CalAsymCol();

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
