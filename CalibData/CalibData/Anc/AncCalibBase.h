// $Header$

/// @file AcdCalibBase
/// @author J. Bogart
#ifndef CalibData_AncCalibBase_h
#define CalibData_AncCalibBase_h

#include <vector>
#include "CalibData/CalibBase.h"
// #include "idents/AcdId.h"

namespace CalibData {
  class AncFinder;
  class RangeBase;

  /**
       Base class for Anc calibration data.  Anc channel
       for tagger is identified by module, layer, chan.
       For qdc we have only module, chan, so internally
       set # layers to 1 and, for generic underpinnings
       expecting a layer number, use zero.

       This class keeps a (pointer to a) vector of pointers to the
       individual per-range datasets and a pointer to a helper class,
       AncFinder, which knows how to compute indices.
  */
  class AncCalibBase : public CalibBase {

  public:
    AncCalibBase(unsigned nModule=1, unsigned nLayer=1, unsigned nChan=32);

    virtual ~AncCalibBase();

    /** 
        Pick out calibration data associated with a particular channel
     */
    virtual RangeBase* getChan(unsigned iMod=0, unsigned iLayer=0, 
                               unsigned iChan=0);

    bool putChan(unsigned iMod, unsigned iLayer, unsigned iChan,
                 RangeBase* data);

    virtual const CLID& clID() const = 0;     // must be overridden
    static const CLID& classID();   // shouldn't get called

    // Maybe won't need to be virtual after all
    virtual StatusCode update(CalibBase& other, MsgStream* log);
    
  protected:
    AncFinder* m_finder;

    std::vector<RangeBase* > m_chans;  // tiles and ribbons

  private:
    static const CLID noCLID;

    /** Due to bug in gcc, gdb can't find symbols in constructors.  This
        method is called by the constructor and does most of the work
    */
    void cGuts(unsigned nMod, unsigned nLay, unsigned nChan);

  };

}  
#endif
