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
       Base class for ACD calibration data, at least for those
       types with a fixed amount of data per tile-pmt-range

       This class keeps a (pointer to a) vector of pointers to the
       individual per-range datasets and a pointer to a helper class,
       AcdFinder, which knows how to compute indices.
  */
  class AcdCalibBase : public CalibBase {

  public:
    AcdCalibBase(const AcdCalibDescription& desc,
		 unsigned nFace=5, unsigned nRow=5, unsigned nCol=5, 
                 unsigned nNA=11, unsigned nPmt=2);

    virtual ~AcdCalibBase();

    virtual const CLID& clID() const = 0;     // must be overridden
    static const CLID& classID();   // shouldn't get called

    // Maybe won't need to be virtual after all
    StatusCode update(CalibBase& other, MsgStream* log);
    
  protected:

    /** 
        Pick out calibration data associated with a particular tile, pmt
     */
    AcdCalibObj* get(idents::AcdId id,unsigned pmt);    
    bool put(idents::AcdId id,unsigned pmt, AcdCalibObj& pmtCalib);
    virtual AcdCalibObj* makeNew() const = 0;

    const AcdCalibDescription* desc() const { return m_desc; }

  private:

    static const CLID noCLID;

    AcdFinder* m_finder;

    std::vector<AcdCalibObj* > m_pmts;  // tiles and ribbons
    std::vector<AcdCalibObj* > m_NAs;   // no detector

    const AcdCalibDescription* m_desc;

    /** Due to bug in gcc, gdb can't find symbols in constructors.  This
        method is called by the constructor and does most of the work
    */
    void cGuts(unsigned nFace, unsigned nRow, unsigned nCol, 
               unsigned nNA, unsigned nPmt);

  };

}  




#endif
