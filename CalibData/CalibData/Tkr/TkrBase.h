// $Header$
/// @file TkrBase.h
/// @author J. Bogart

#ifndef CalibData_TkrBase_h
#define CalibData_TkrBase_h

#include <vector>
#include "CalibData/CalibBase.h"
#include "CalibData/Tkr/TkrFinder.h"
#include "idents/TkrId.h"

namespace CalibData {
  /**
     Base class for Tkr calibrations other than bad strips, which have
     a fixed amount of data per Si layer or fe chip.
  */
  class RangeBase;
  class TkrBase : public CalibBase {
    
  public:
    /**
       Constructor configures its TkrFinder and keeps track of whether
       data is stored directly in a vector or whether its a vector
       of pointers.
    */
    TkrBase(unsigned nTowerRow=4, unsigned nTowerCol=4, unsigned nTray=19,
            unsigned nFeChip=0, bool indirect=true);
    virtual ~TkrBase();
    
    virtual RangeBase* getChannel(const idents::TkrId& id, unsigned feChip=0);

    virtual RangeBase* getChannel(unsigned towerRow, unsigned towerCol,
                                  unsigned tray, bool top, unsigned feChip=0);
    
    virtual bool putChannel(RangeBase* data, const idents::TkrId& id, 
                            unsigned feChip=0);

    virtual bool putChannel(RangeBase* data, unsigned towerRow, 
                            unsigned towerCol, unsigned tray, 
                            bool top, unsigned feChip=0);

    // Get dimensioning information; needed when transforming to 
    // permanent storage
    /// Get # tower rows
    unsigned getNTowerRow() const {return m_finder->getNTowerRow();}

    /// Get # tower columns
    unsigned getNTowerCol() const {return m_finder->getNTowerCol();}

    /// Get #  trays
    unsigned getNUnilayer() const {return m_finder->getNUnilayer();}

    /// Get # fe chips / unilayer
    unsigned getNChip() const {return m_finder->getNChip();}
    virtual const CLID& clID() const = 0;     // must be overridden
    static const CLID& classID();   // shouldn't get called

    // Maybe won't need to be virtual after all
    virtual StatusCode update(CalibBase& other, MsgStream* log);
    
  protected:
    TkrFinder* m_finder;
    //    std::vector<RangeBase* >* m_pR;
    std::vector<RangeBase* > m_ranges;

    /// Default is true: keep vector of pointers to data.  Else derived
    /// class keeps vector of data values and must do its own fetching
    /// and putting.
    bool m_indirect;
    
    // cache last index found, for use of derived classes
    unsigned m_ix;
    bool     m_ixValid;

  private:
    static const CLID noCLID;

    /** Due to bug in gcc, gdb can't find symbols in constructors.  This
        method is called by the constructor and does most of the work
    */
    void cGuts(unsigned nTowerRow, unsigned nTowerCol, 
               unsigned nTray, unsigned nChip);
  };

}  
#endif

