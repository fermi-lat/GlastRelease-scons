// $Header$

/// @file CalTholdMuon.h
/** 
  This file contains class definitions corresponding to the calibration
  type  CAL_TholdMuon.    Two classes are defined.  One, CalTholdMuon, 
  is the data pertaining to a  single crystal face. The other, 
  CalTholdMuonCol, is a collection of CalTholdMuon instances. 

 @author J. Bogart
*/

#ifndef CalibData_CalTholdMuon_h
#define CalibData_CalTholdMuon_h

#include "CalibData/RangeBase.h"
#include "CalibData/Cal/CalCalibBase.h"

namespace CalibData {
  /** @class CalTholdMuon
     Data pertaining to a single crystal face.  It inherits from RangeBase.
  */
  class CalTholdMuon : public RangeBase {
  public:
    /** 
      The ped vector is indexed using constants defined in the enum
      idents::CalXtalId::AdcRange. 
     */
    CalTholdMuon(const ValSig* FLE = 0, 
                 const ValSig* FHE = 0,
                 const std::vector<ValSig>* ped = 0);
    ~CalTholdMuon();

    // Bunch of get functions

    const ValSig* getFLE() const {return &m_FLE;}
    const ValSig* getFHE() const {return &m_FHE;}

    /**
       Returns pedestal, uncertainty for a single range
       Caller should use enum definitions in idents::CalXtalId::AdcRange
       for @a range argument. If @a range is not one of these enum
       values or if specified range has no data, method will return
       a null pointer.
    */
    const ValSig* getPed(int range) const;

    /// Returns vector of peds for all ranges.
    const std::vector<ValSig>* getPeds() const {return m_ped; }

    virtual void update(RangeBase* other);

  private:
    /// Threshold for fast low energy shaper for muon data
    ValSig m_FLE;
    /// Threshold for fast high energy shaper for muon data
    ValSig m_FHE;

    /// Pedestals as measured during muon calibration, indexed a la
    /// idents::CalXtalId::AdcRange
    std::vector<ValSig>* m_ped;
  };

  /**  @class CalTholdMuonCol
     This class represents CalTholdMuon information for the entire detector.
     Individual CalTholdMuon instances may be looked up by CalXtalId.
  */
  class CalTholdMuonCol :  public CalCalibBase {

  public:
    CalTholdMuonCol(unsigned nTowerRow=4, unsigned nTowerCol=4, 
                    unsigned nLayer=8, unsigned nXtal=12, unsigned nFace=2);

    ~CalTholdMuonCol();

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
