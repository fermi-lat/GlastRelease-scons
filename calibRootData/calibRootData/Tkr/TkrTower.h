#ifndef calibRootData_TkrTower_h
#define calibRootData_TkrTower_h

/**
 *  @file TkrTower.h
 *  @brief   trivial class to identify tracker tower
 */

#include "TObject.h"
#include "TString.h"

namespace calibRootData {
  /** 
   *    @class TkrTower.h
   *    @brief Identify by location and/or serial number
   */
  class TkrTower : public TObject {
  private:
    UInt_t m_row;
    UInt_t m_col;
    TString m_hwserial;

  public:
    TkrTower(UInt_t row=0, UInt_t col=0) : m_row(row), m_col(col) {
      m_hwserial = TString();
    }
    TkrTower(UInt_t row, UInt_t col, TString hwserial) : 
      m_row(row), m_col(col) {
      m_hwserial = hwserial;
    }
    ~TkrTower() {}

    void init(UInt_t row, UInt_t col, TString hwserial) {
      m_row = row;
      m_col = col;
      m_hwserial = hwserial;
    }
    UInt_t getRow() {return m_row;}
    UInt_t getCol() {return m_col;}
    TString getSerial() {return m_hwserial;}
    
  private:    

    ClassDef(calibRootData::TkrTower, 1)

 };
}

#endif
