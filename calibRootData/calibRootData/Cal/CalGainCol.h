#ifndef calibRootData_CalGainCol_h
#define calibRootData_CalGainCol_h

#include "TObject.h"
#include "TString.h"
#include "../Cal/CalBase.h"
#include "../Cal/CalGain.h"
/*  Might need some more ROOT includes
#include "TRef.h"
#include "TRefArray.h"
*/

/** @class calibRootData::CalBase
* @brief Base class for (almost) all Cal calibration ROOT classes.  
* 
* Those representing hot/dead channels probably won't use this base.
*
* @author Joanne Bogart
*  
* $Header$
*/
namespace calibRootData {
  class CalGainCol : public CalBase  {
    
  public:


    CalGainCol(TString instr = TString(), TString t = TString(),
               unsigned nRow=0, unsigned nCol=0, unsigned nLayer=0, 
               unsigned nXtal=0, unsigned nFace=0, unsigned nRange=0) : 
      CalBase(i, t, nRow, nCol, nLayer, nXtal, nFace, nRange)
       {m_gainObj.clear() };
    virtual ~CalGainCol() {m_gainObj.clear };

    // Will this work as it stands or do I need to make a (non-const) copy?
    void  addChannel(const CalGain& gainObj) {gains.push_back(gainObj);}

    void addChannel(const commonRootData::CalXtalId id, Float_t gain, 
                    Float_t sig) {
      gains.push_back(CalGain(id.getPackedId(), gain, sig));
    }
    
  private:

    std::vector<CalGain> gains;

    ClassDef(calibRootData::CalGainCol, 1)
  };
}
#endif
