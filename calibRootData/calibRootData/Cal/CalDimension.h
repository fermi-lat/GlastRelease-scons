#ifndef calibRootData_CalDimension_h
#define calibRootData_CalDimension_h

#include "TObject.h"
#include "TString.h"
/*  Might need some more ROOT includes
#include "TRef.h"
#include "TRefArray.h"
*/

/** @class calibRootData::CalDimension
* @brief Base class for (almost) all Cal calibration ROOT classes.  
* 
* Those representing hot/dead channels probably won't use this base.
*
* @author Joanne Bogart
*  
* $Header$
*/
namespace calibRootData {
  class CalDimension : public TObject {
  public:
    
    CalDimension(unsigned nRow=0, unsigned nCol=0, 
                 unsigned nLayer=0, unsigned nXtal=0, unsigned nFace=0, 
                 unsigned nRange=0, unsigned nDacCol=0, bool exact=false) : 
      m_nRow(nRow), m_nCol(nCol), m_nLayer(nLayer),
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

 protected:

    UChar_t m_nRow;
    UChar_t m_nCol;
    UChar_t m_nLayer;
    UChar_t m_nXtal;
    UChar_t m_nFace;
    UChar_t m_nRange;
    UChar_t m_nDacCol;
    Bool_t  m_exact;    

  private:
    friend class CalBase;  // so it can call initialize

    void initialize(unsigned nRow, unsigned nCol, unsigned nLayer, 
                    unsigned nXtal, unsigned nFace, unsigned nRange, 
                    unsigned nDacCol) {
      m_nRow = nRow;    m_nCol = nCol;  m_nLayer = nLayer;   m_nXtal = nXtal;
      m_nFace = nFace;  m_nRange = nRange;      m_nDacCol = nDacCol;
    }
  
    ClassDef(calibRootData::CalDimension, 1)
  };
}
#endif
