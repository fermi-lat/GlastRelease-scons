// $Header$

/// @file CalTholdCI.h
/** 
  This file contains class definitions corresponding to the calibration
  type  CAL_TholdCI.    Two classes are defined.  One, CalTholdCI, 
  is the data pertaining to a  single crystal face. The other, CalTholdCICol, 
  manages the full collection of CalTholdCI instances for the detector

 @author J. Bogart
*/

#ifndef CalibData_CalTholdCI_h
#define CalibData_CalTholdCI_h

#include "CalibData/RangeBase.h"
#include "CalibData/Cal/CalCalibBase.h"

namespace CalibData {
  /** @class CalTholdCI
     Data pertaining to a single crystal face.  It inherits from RangeBase.
  */
  class CalTholdCI : public RangeBase {
  public:
    CalTholdCI(const std::vector<ValSig>* ULD = 0,
               const ValSig* FLE = 0, 
               const ValSig* FHE = 0,
               const ValSig* LAC = 0,
               const std::vector<ValSig>* ped = 0);
    ~CalTholdCI();


    /**
       Returns ULD, uncertainty for a single range.
       Caller should use enum definition in idents::CalXtalId::AdcRange
       for @a range argument.  If @a range is not one of these enum
       values or if specified range has no data, method will return
       a null pointer.
     */
    const ValSig* getULD(int range) const;

    /// returns ULD threshold for all ranges
    const std::vector<ValSig>* getULDs() const {return m_ULD;}

    const ValSig* getFLE() const {return &m_FLE;}
    const ValSig* getFHE() const {return &m_FHE;}
    const ValSig* getLAC() const {return &m_LAC;}

    /**
       Returns pedestal, uncertainty for a single range
       Caller should use enum definition in idents::CalXtalId::AdcRange
       for @a range argument.  If @a range is not one of these enum
       values or if specified range has no data, method will return
       a null pointer.
     */
    const ValSig* getPed(int range) const;

    // Returns a vector of peds for all ranges
    const std::vector<ValSig>* getPeds() const {return m_ped;}

    virtual void update(RangeBase* other);

  private:
    ///  Upper level discriminator per-range threshold, indexed by
    ///  idents::CalXtalId::AdcRange
    std::vector<ValSig>* m_ULD;
    /// threshold for the fast low energy shaper for charge injection data
    ValSig m_FLE;
    /// threshold for the fast high energy shaper for charge injection data
    ValSig m_FHE;
    /// threshold for the log accept discriminator for charge injection data
    ValSig m_LAC;
    /// Pedestals as measured during charge injection calibration, indexed
    /// a la idents::CalXtalId::AdcRange
    std::vector<ValSig>* m_ped;
  };

  /**  @class CalTholdCICol
     This class represents CalTholdCI information for the entire detector.
     Individual CalTholdCI instances may be looked up by CalXtalId.
  */
  class CalTholdCICol :  public CalCalibBase {

  public:
    CalTholdCICol(unsigned nTowerRow=4, unsigned nTowerCol=4, 
                    unsigned nLayer=8, unsigned nXtal=12, unsigned nFace=2);

    ~CalTholdCICol();

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
