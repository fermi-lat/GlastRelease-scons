#ifndef calibRootData_CalDimension_h
#define calibRootData_CalDimension_h

#include "TObject.h"
#include "TString.h"

/** @class calibRootData::CalDimension
* @brief Dimension information needed by TDS to index Calorimeter
* calibration data.  Also includes instrument name.
*
* @author Joanne Bogart
*  
* $Header$
*/
namespace calibRootData {
  class CalDimension : public TObject {
  public:
    
    CalDimension(TString instr= TString(), unsigned nRow=0, unsigned nCol=0, 
                 unsigned nLayer=0, unsigned nXtal=0, unsigned nFace=0, 
                 unsigned nRange=0, unsigned nDacCol=0, bool exact=false) : 
      m_instrument(instr), m_nRow(nRow), m_nCol(nCol), m_nLayer(nLayer),
      m_nXtal(nXtal), m_nFace(nFace), m_nRange(nRange), m_nDacCol(nDacCol),
      m_exact(exact) { };

    virtual ~CalDimension() { };

    unsigned getNRow() const {return m_nRow;}
    unsigned getNCol() const {return m_nCol;}
    unsigned getNLayer() const {return m_nLayer;}
    unsigned getXtal() const {return m_nXtal;}
    unsigned getNFace() const {return m_nFace;}
    unsigned getNRange() const {return m_nRange;}
    unsigned getNDacCol() const {return m_nDacCol;}
    bool     isExact() const {return m_exact;}

    /// return true if arg. matches internal instrument name
    bool isInstrument(const TString instrument) const {
      return (instrument == m_instrument);
    }

    /// Set all dimension variables (but not instrument name)
    void initialize(unsigned nRow, unsigned nCol, unsigned nLayer, 
                    unsigned nXtal, unsigned nFace, unsigned nRange, 
                    unsigned nDacCol) {
      m_nRow = nRow;    m_nCol = nCol;  m_nLayer = nLayer;   m_nXtal = nXtal;
      m_nFace = nFace;  m_nRange = nRange;      m_nDacCol = nDacCol;
    }
    void setInstrument(TString i) {m_instrument = i;}
 private:

    TString m_instrument;
    UChar_t m_nRow;
    UChar_t m_nCol;
    UChar_t m_nLayer;
    UChar_t m_nXtal;
    UChar_t m_nFace;
    UChar_t m_nRange;
    UChar_t m_nDacCol;
    Bool_t  m_exact;    

    ClassDef(calibRootData::CalDimension, 1)
  };
}
#endif
