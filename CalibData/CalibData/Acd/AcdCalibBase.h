// $Header$

/// @file AcdCalibBase
/// @author J. Bogart
#ifndef CalibData_AcdCalibBase_h
#define CalibData_AcdCalibBase_h

#include <vector>
#include "CalibData/CalibBase.h"
#include "CalibData/Acd/AcdCalibEnum.h"
#include "idents/AcdId.h"

namespace CalibData {
  class AcdFinder;
  class AcdCalibObj;
  class AcdCalibDescription;

  /**
   * @class AcdCalibBase
   * 
   * @brief Base class for ACD calibration data sets.
   * 
   * This class keeps a (pointer to a) vector of pointers to the
   * individual channel calibrations and a pointer to a helper class,
   * AcdFinder, which knows how to compute indices.
   *
   **/

  class AcdCalibBase : public CalibBase {

  public:
    /// Build space for the whole ACD
    AcdCalibBase(const AcdCalibDescription& desc,
		 unsigned nFace=5, unsigned nRow=5, unsigned nCol=5, 
                 unsigned nNA=11, unsigned nPmt=2);

    /// Cleanup
    virtual ~AcdCalibBase();

    virtual const CLID& clID() const = 0;   // must be overridden  
    static const CLID& classID();           // shouldn't get called

    /// Update all the values
    StatusCode update(CalibBase& other, MsgStream* log);
    
  protected:

    /// Pick out calibration data associated with a particular tile, pmt
    AcdCalibObj* get(idents::AcdId id,unsigned pmt);    

    /// Put in  alibration data associated with a particular tile, pmt
    bool put(idents::AcdId id,unsigned pmt, AcdCalibObj& pmtCalib);

    /// Allocate space for calibration data for one channel
    virtual AcdCalibObj* makeNew() const = 0;

    /// return the descrption of the calibration this object manages
    const AcdCalibDescription* desc() const { return m_desc; }

  private:

    static const CLID noCLID;

    /// This guy knows how to compute indices.
    AcdFinder* m_finder;

    /// tiles and ribbons
    std::vector<AcdCalibObj* > m_pmts;
    /// unattached channels
    std::vector<AcdCalibObj* > m_NAs;

    /// the descrption of the calibration this object manages
    const AcdCalibDescription* m_desc;

    /** Due to bug in gcc, gdb can't find symbols in constructors.  This
        method is called by the constructor and does most of the work
    */
    void cGuts(unsigned nFace, unsigned nRow, unsigned nCol, 
               unsigned nNA, unsigned nPmt);

  };

}  




#endif
