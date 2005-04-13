#ifndef calibRootData_CalGainCol_h
#define calibRootData_CalGainCol_h

#include "../Cal/CalGain.h"
#include "../Cal/CalDimension.h"

/** @class calibRootData::CalGainCol
* @brief Calorimeter gains calibration
*
* @author Joanne Bogart
*  
* $Header$
*/
namespace calibRootData {
  class CalGainCol : public TObject  {
    
  public:


    CalGainCol(){}
    virtual ~CalGainCol() {m_gains.clear(); };

    // Will this work as it stands or do I need to make a (non-const) copy?
    void  addChannel(const CalGain& gainObj) {m_gains.push_back(gainObj);}

    void addChannel(const commonRootData::CalXtalId id, Float_t gain, 
                    Float_t sig) {
      m_gains.push_back(CalGain(id.getPackedId(), gain, sig));
    }

    ///  Transparent access to CalDimension part of class.  
    CalDimension* getDimension() {return &m_dim;}
    std::vector<calibRootData::CalGain>& getGains() {return m_gains;}
    
  private:

    std::vector<CalGain> m_gains;
    CalDimension         m_dim;
   
    ClassDef(calibRootData::CalGainCol, 1)
  };
}
#endif
