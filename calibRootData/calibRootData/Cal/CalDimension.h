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
                 unsigned nLayer=0, unsigned nXtal=0, 
                 unsigned nFace=0, unsigned nRange=0) : 
      m_nRow(nRow), m_nCol(nCol), m_nLayer(nLayer),
      m_nXtal(nXtal), m_nFace(nFace), m_nRange(nRange) { };

    virtual ~CalDimension() { };

    unsigned getNRow() const {return m_nRow;}
    unsigned getNCol() const {return m_nCol;}
    unsigned getNLayer() const {return m_nLayer;}
    unsigned getXtal() const {return m_nXtal;}
    unsigned getNFace() const {return m_nFace;}
    unsigned getNRange() const {return m_nRange;}

 protected:

    unsigned char m_nRow;
    unsigned char m_nCol;
    unsigned char m_nLayer;
    unsigned char m_nXtal;
    unsigned char m_nFace;
    unsigned char m_nRange;
  
    ClassDef(calibRootData::CalDimension, 1)
  };
}
