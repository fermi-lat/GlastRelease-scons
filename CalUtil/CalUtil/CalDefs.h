#ifndef CalDefs_H
#define CalDefs_H
// $Header$

// LOCAL

// GLAST
#include "idents/CalXtalId.h"

// EXTLIB

// STD
#include <string>
#include <vector>
#include <ostream>
#include <stdexcept>
#include <sstream>

using namespace std;

/** @file CalDefs.h
    \author Zachary Fewtrell

    \brief enums & array index types for use throughout GLAST Cal software.  
    Required mainly because idents::CalXtalId is non-contiguous bitfield structure which
    can't (reasonably) be used as an array index.

    All classes have the following properties:

    - Internal integer representation is contiguous, so they are suitable for array
    indexing.
    
    - Iteration operators defined for use in loops.  
    
    - Full set of conversion routines defined for converting between related indices.
    
    - All classes are wrappers for simple integers, can efficiently be placed on function stack.
    
    - Classes exist for indexing components within a single crystal or througout
    the entire LAT.  Appropriate conversion routines are provided.

*/
namespace CalUtil {
  // forward declarations (sometimes needed)
  class DirNum;
  class RngNum;

    
  /// generic type2string converter
  template <typename _T>
  std::string toString(const _T &val) {
    std::ostringstream tmp;
    tmp << val;
    return tmp.str();
  }

  /** \brief Atomic (simple) Cal Geometry Id's

      Id single Cal component w/ in it's immediate container
  */
  class SimpleId {
  public:
    /// prefix ++ operator
    SimpleId operator++() {
      SimpleId tmp(*this);
      m_data++;
      return tmp;
    }

    /// postfix ++ operator
    SimpleId& operator++(const int) {
      m_data++; 
      return *this;
    }

    unsigned short val() const {return m_data;}

    std::string toStr() const {return toString(m_data);}

  protected:
    explicit SimpleId(const unsigned short val=0) : m_data(val) {}
    unsigned short m_data;
  };

  /// index class for GLAST Cal module tower bay
  class TwrNum : public SimpleId {
  public:
    /// \note not explicit b/c 0-15 is unabiguous, well known index
    TwrNum(const unsigned short val=0) : SimpleId(val) {}

    TwrNum(const unsigned short tRow, const unsigned short tCol) : 
      SimpleId(tRow*N_COLS + tCol) {}

    unsigned short  getRow() const {return m_data/N_COLS;}
    unsigned short  getCol() const {return m_data%N_COLS;}

    static const unsigned short N_VALS = 16; 
    bool isValid() const {return m_data < N_VALS;}
            
    static const unsigned short N_COLS = 4;
    static const unsigned short N_ROWS = 4;

    bool operator==(const TwrNum that) const {return m_data == that.m_data;}
    bool operator!=(const TwrNum that) const {return m_data != that.m_data;}
    bool operator>=(const TwrNum that) const {return m_data >= that.m_data;}
    bool operator<=(const TwrNum that) const {return m_data <= that.m_data;}
    bool operator>(const TwrNum that) const {return m_data > that.m_data;}
    bool operator<(const TwrNum that) const {return m_data < that.m_data;}

  };

  /// online TEM indexing is identical to offline Tower indexing.
  typedef TwrNum GTEMNum;

  class LyrNum;

  /// index class for GLAST Cal GCRC (layer 0-3 on a single cal front end board).
  class GCRCNum : public SimpleId {
  public:
    explicit GCRCNum(const unsigned short val=0) : SimpleId(val) {}

    GCRCNum(const LyrNum &lyr);

    static const unsigned short N_VALS=4;
    bool isValid() const {return m_data < N_VALS;}

    bool operator==(const GCRCNum that) const {return m_data == that.m_data;}
    bool operator!=(const GCRCNum that) const {return m_data != that.m_data;}
    bool operator>=(const GCRCNum that) const {return m_data >= that.m_data;}
    bool operator<=(const GCRCNum that) const {return m_data <= that.m_data;}
    bool operator>(const GCRCNum that) const {return m_data > that.m_data;}
    bool operator<(const GCRCNum that) const {return m_data < that.m_data;}
  };

  /// index class for GLAST Cal xtal layer w/in Cal module
  class LyrNum : public SimpleId {
  public:
    /// \note explicit b/c of alternate layer numbering schemes
    explicit LyrNum(const unsigned short val=0) : SimpleId(val) {}
    LyrNum(const DirNum dir, const GCRCNum gcrc);
    
    static const unsigned short N_VALS=8;
    bool isValid() const {return m_data < N_VALS;}

    // 0 is 'X' Direction
    DirNum getDir() const;
    GCRCNum getGCRC() const {return GCRCNum(m_data/2);}

    bool operator==(const LyrNum that) const {return m_data == that.m_data;}
    bool operator!=(const LyrNum that) const {return m_data != that.m_data;}
    bool operator>=(const LyrNum that) const {return m_data >= that.m_data;}
    bool operator<=(const LyrNum that) const {return m_data <= that.m_data;}
    bool operator>(const LyrNum that) const {return m_data > that.m_data;}
    bool operator<(const LyrNum that) const {return m_data < that.m_data;}

  };


  /// index class for GLAST Cal xtal direction ('X' or 'Y')
  class DirNum : public SimpleId {
  public:
    explicit DirNum(const unsigned short val=0) : SimpleId(val) {}

    static const unsigned short N_VALS=2;
    bool isValid() const {return m_data < N_VALS;}

    bool operator==(const DirNum that) const {return m_data == that.m_data;}
    bool operator!=(const DirNum that) const {return m_data != that.m_data;}

    const std::string &toStr() const;

  };
  const DirNum X_DIR(0);
  const DirNum Y_DIR(1);

  /// index class for GLAST Cal xtal column w/in Cal xtal layer
  class ColNum : public SimpleId {
  public:
    /// \note not explicit b/c 0-11 unambiguous index
    ColNum(const unsigned short val=0) : SimpleId(val) {}

    static const unsigned short N_VALS=12;
    bool isValid() const {return m_data < N_VALS;}

    bool operator==(const ColNum that) const {return m_data == that.m_data;}
    bool operator!=(const ColNum that) const {return m_data != that.m_data;}
    bool operator>=(const ColNum that) const {return m_data >= that.m_data;}
    bool operator<=(const ColNum that) const {return m_data <= that.m_data;}
    bool operator>(const ColNum that) const {return m_data > that.m_data;}
    bool operator<(const ColNum that) const {return m_data < that.m_data;}

  };

  /// GCFE numbering is identical to offline 'column' indexing
  typedef ColNum GCFENum;

  /// index class for GLAST Cal xtal face ('POS' or 'NEG')
  class FaceNum : public SimpleId {
  public:
    /// \note explicit b/c face indexing is ambiguous 
    explicit FaceNum(const idents::CalXtalId::XtalFace face) : SimpleId(face) {}

    /// \note explicit b/c face indexing is ambiguous
    explicit FaceNum(const unsigned short val=0) : SimpleId(val) {}

    /// inverse of toStr() method
    FaceNum(const std::string &str);
    
    static const vector<string> MNEM;
    
    static const unsigned short N_VALS=2;    
    
    bool isValid() const {return m_data < N_VALS;}
    
    const std::string &toStr() const;
    
    /// allow quick conversion to idents::CalXtalId::XtalFace since internal
    /// storage is the same
    operator idents::CalXtalId::XtalFace() const {
      return (idents::CalXtalId::XtalFace)m_data;
    }

    bool operator==(const FaceNum that) const {return m_data == that.m_data;}
    bool operator!=(const FaceNum that) const {return m_data != that.m_data;}

    /// return face id on opposite side of xtal
    FaceNum oppositeFace() const {
      return FaceNum(1-m_data);
    }
  };

  const FaceNum POS_FACE(idents::CalXtalId::POS);
  const FaceNum NEG_FACE(idents::CalXtalId::NEG);

  /// index class for GLAST Cal GCCC (face 0-3 on a single cal TEM module).
  /// 0=X+,1=Y+,2=X-,3=Y-
  class GCCCNum : public SimpleId {
  public:
    explicit GCCCNum(const unsigned short val=0) : SimpleId(val) {}

    GCCCNum(const DirNum dir,
            const FaceNum face) 
    {
      m_data = 0;
      if (face == NEG_FACE)
        m_data += 2;
      if (dir == Y_DIR)
        m_data++;
    }
      

    static const unsigned short N_VALS=4;
    bool isValid() const {return m_data < N_VALS;}

    bool operator==(const GCCCNum that) const {return m_data == that.m_data;}
    bool operator!=(const GCCCNum that) const {return m_data != that.m_data;}
    bool operator>=(const GCCCNum that) const {return m_data >= that.m_data;}
    bool operator<=(const GCCCNum that) const {return m_data <= that.m_data;}
    bool operator>(const GCCCNum that) const {return m_data > that.m_data;}
    bool operator<(const GCCCNum that) const {return m_data < that.m_data;}
  };

  const GCCCNum X_POS(X_DIR, POS_FACE);
  const GCCCNum Y_POS(Y_DIR, POS_FACE);
  const GCCCNum X_NEG(X_DIR, NEG_FACE);
  const GCCCNum Y_NEG(Y_DIR, NEG_FACE);



  /// index class for GLAST Cal photo-diode on single xtal face 
  /// ('LARGE' or 'SMALL')
  class DiodeNum : public SimpleId {
  public:
    explicit DiodeNum(const unsigned short val=0) : SimpleId(val) {}

    /// inverse of toStr() method
    DiodeNum(const std::string &str);

    RngNum getX8Rng() const;
    RngNum getX1Rng() const;

    static const vector<string> MNEM;
    static const unsigned short N_VALS = 2; 
    bool isValid() const {return m_data < N_VALS;}

    const std::string &toStr() const;

    friend ostream &operator<<(ostream &stream, const DiodeNum &id);

    bool operator==(const DiodeNum that) const {return m_data == that.m_data;}
    bool operator!=(const DiodeNum that) const {return m_data != that.m_data;}


  };
  const DiodeNum LRG_DIODE(idents::CalXtalId::LARGE);
  const DiodeNum SM_DIODE (idents::CalXtalId::SMALL);

  ostream &operator<<(ostream &stream, const DiodeNum &id);

  /// THX (Track & Hold Multiplier) can be either X8 or X1
  class THXNum : public SimpleId {
  public:
    explicit THXNum(const unsigned short val=0) : SimpleId(val) {}

    static const vector<string> MNEM;

    bool isValid() const {return m_data < N_VALS;}
    static const unsigned short N_VALS=2;
    const std::string &toStr() const;

    bool operator==(const THXNum that) const {return m_data == that.m_data;}
    bool operator!=(const THXNum that) const {return m_data != that.m_data;}

  };
  const THXNum THX8(0);
  const THXNum THX1(1);

  /// index class for GLAST Cal ADC range (LEX8 -> HEX1)
  class RngNum : public SimpleId {
  public:
    /// \not not explicit b/c 0-3 is unambiguous range indexing scheme
    RngNum(const unsigned short val=0) : SimpleId(val) {}

    RngNum(DiodeNum diode, THXNum thx) :
      SimpleId(diode.val()*2 + thx.val()) {}

    static const vector<string> MNEM;
    static const unsigned short N_VALS=4;
    const std::string &toStr() const;
    
    DiodeNum getDiode() const {
      using idents::CalXtalId;
      return DiodeNum(idents::CalXtalId::rangeToDiode((idents::CalXtalId::AdcRange)m_data));
    }
    bool isValid() const {return m_data < N_VALS;}

    operator idents::CalXtalId::AdcRange() const {
      return (idents::CalXtalId::AdcRange)m_data;
    }

    THXNum getTHX() const {return THXNum(m_data%2);}

    bool operator==(const RngNum that) const {return m_data == that.m_data;}
    bool operator!=(const RngNum that) const {return m_data != that.m_data;}
    bool operator>=(const RngNum that) const {return m_data >= that.m_data;}
    bool operator<=(const RngNum that) const {return m_data <= that.m_data;}
    bool operator>(const RngNum that) const {return m_data > that.m_data;}
    bool operator<(const RngNum that) const {return m_data < that.m_data;}

  };

  const RngNum LEX8(idents::CalXtalId::LEX8);
  const RngNum LEX1(idents::CalXtalId::LEX1);
  const RngNum HEX8(idents::CalXtalId::HEX8);
  const RngNum HEX1(idents::CalXtalId::HEX1);


  /** \brief Id all components of same type within a single CAL xtal.
  
      internal integer storage is contiguous and can be used as an array Index
  */

  class XtalWideIndex {
  public:
    /// prefix ++ operator
    XtalWideIndex operator++() {
      XtalWideIndex tmp(*this);
      m_data++;
      return tmp;
    }

    /// postfix ++ operator
    XtalWideIndex& operator++(const int) {
      m_data++; 
      return *this;
    }

    bool operator==(const XtalWideIndex that) const {
      return m_data == that.m_data;}
    bool operator<=(const XtalWideIndex that) const {
      return m_data <= that.m_data;}
    bool operator!=(const XtalWideIndex that) const {
      return m_data != that.m_data;}
    bool operator< (const XtalWideIndex that) const {
      return m_data <  that.m_data;}

    unsigned val() const {return m_data;}

  protected:
    explicit XtalWideIndex(const unsigned short val=0) : m_data(val) {}
    unsigned short m_data;
  };

  /// index class for all 4 photo diodes in one Cal xtal.
  class XtalDiode : public XtalWideIndex {
  public: 
    XtalDiode(const FaceNum face, const DiodeNum diode) :
      XtalWideIndex(face.val()*FACE_BASE + diode.val()) {}
    XtalDiode() : XtalWideIndex() {}

    DiodeNum getDiode() const {return DiodeNum(m_data%FACE_BASE);}
    FaceNum getFace()  const {return FaceNum((idents::CalXtalId::XtalFace)(m_data/FACE_BASE));}

    static const unsigned short N_VALS = FaceNum::N_VALS*DiodeNum::N_VALS;
    bool isValid() const {return m_data < N_VALS;}
  protected:
    static const unsigned short FACE_BASE = DiodeNum::N_VALS;
  };

  /// index class for all 8 ADC readouts on one Cal xtal.
  class XtalRng : public XtalWideIndex {
  public:
    XtalRng(const FaceNum face, const RngNum rng) :
      XtalWideIndex(face.val()*FACE_BASE + rng.val()) {};
    XtalRng() : XtalWideIndex() {}

    FaceNum getFace() const {return FaceNum(idents::CalXtalId::XtalFace(m_data/FACE_BASE));}
    RngNum  getRng() const {return RngNum(m_data%FACE_BASE);}

    XtalDiode getXtalDiode() const {
      return XtalDiode(getFace(), getRng().getDiode());
    }

    static const unsigned short N_VALS = FaceNum::N_VALS*RngNum::N_VALS;
    bool isValid() const {return m_data < N_VALS;}
  protected:
    static const unsigned short FACE_BASE = RngNum::N_VALS;
  };

  /////////////////////////////////////////////
  /// Tower Wide Cal Indexes                ///
  /// Id all components of same type within ///
  /// single tower.
  /// internal integer storage is contiguous///
  /// and can be used as an array Index     ///
  /////////////////////////////////////////////

  /** \brief Abstract root class for indexing Cal 
      components w/in single tower
    
      any interfaces w/ idents::CalXtalId always work w/ tower bay 
      # 0.
  */
  class TwrWideIndex {
  public:  
    /// prefix ++ operator
    TwrWideIndex operator++() {
      TwrWideIndex tmp(*this);
      m_data++;
      return tmp;
    }
   
    /// postfix ++ operator
    TwrWideIndex& operator++(int) {
      m_data++; 
      return *this;
    }
     
    unsigned val() const {return m_data;}
       
    bool operator==(const TwrWideIndex that) const {return m_data == that.m_data;}
    bool operator<=(const TwrWideIndex that) const {return m_data <= that.m_data;}
    bool operator!=(const TwrWideIndex that) const {return m_data != that.m_data;}
    bool operator< (const TwrWideIndex that) const {return m_data <  that.m_data;}
    //TwrWideIndex& operator= (const TwrWideIndex that) {m_data = that.m_data;}
  protected:
    explicit TwrWideIndex(const unsigned val=0) : m_data(val) {}
    unsigned m_data;
  };

  /// index all crytals in single cal module
  class tXtalIdx : public TwrWideIndex {
  public:
    explicit tXtalIdx(const idents::CalXtalId xtal) :
      TwrWideIndex(calc(LyrNum(xtal.getLayer()),
                        xtal.getColumn())) {}
   
    tXtalIdx(const LyrNum lyr, const ColNum col) :
      TwrWideIndex(calc(lyr, col)) {}
       
    tXtalIdx() : TwrWideIndex() {}
       
    static const unsigned N_VALS  = LyrNum::N_VALS*ColNum::N_VALS;
           
    LyrNum getLyr() const {return LyrNum((m_data)/LYR_BASE);}
    ColNum getCol() const {return m_data%LYR_BASE;}
               
    bool isValid() const {return m_data < N_VALS;}
                   
  private:
    static unsigned calc(const LyrNum lyr, const ColNum col) {
      return lyr.val()*LYR_BASE + col.val();
    }
    static const unsigned LYR_BASE  = ColNum::N_VALS;
  };

  /// index all crystal faces in single cal module
  class tFaceIdx : public TwrWideIndex {
  public:
    tFaceIdx() : TwrWideIndex() {}
   
    explicit tFaceIdx(const idents::CalXtalId xtalId) {
      if (!xtalId.validFace())
        throw invalid_argument("tFaceIdx requires valid face info in xtalId."
                               "  Programmer error");
                           
      m_data = calc(LyrNum(xtalId.getLayer()),
                    xtalId.getColumn(),
                    FaceNum((idents::CalXtalId::XtalFace)xtalId.getFace()));
    }
     
    tFaceIdx(const LyrNum lyr,const  ColNum col, const FaceNum face) :
      TwrWideIndex(calc(lyr,col,face)) {}
       
    tFaceIdx(const tXtalIdx xtal, const FaceNum face) :
      TwrWideIndex(calc(xtal.getLyr(), xtal.getCol(), face)) {}
         
    static const unsigned N_VALS = tXtalIdx::N_VALS*FaceNum::N_VALS;
             
    tXtalIdx getTXtalIdx() const {return tXtalIdx(getLyr(),
                                                  getCol());}
               
    LyrNum getLyr()  const {return LyrNum((m_data)/LYR_BASE);}
    ColNum getCol()  const {return (m_data%LYR_BASE)/COL_BASE;}
    FaceNum getFace() const {return FaceNum((idents::CalXtalId::XtalFace)(m_data%COL_BASE));}
                                         
    bool isValid() const {return m_data < N_VALS;}
  private:
    static unsigned calc(const LyrNum lyr, const ColNum col, const FaceNum face) {
      return lyr.val()*LYR_BASE + col.val()*COL_BASE + face.val();
    }
    static const unsigned short COL_BASE  = FaceNum::N_VALS;
    static const unsigned short LYR_BASE  = COL_BASE*ColNum::N_VALS;
  };

  /// index all crystal diodes in single cal module
  class tDiodeIdx : public TwrWideIndex {
  public:
    tDiodeIdx(const idents::CalXtalId xtalId, DiodeNum diode) {
      if (!xtalId.validFace())
        throw invalid_argument("tDiodeIdx requires valid face info in xtalId.  Programmer error");
                           
      m_data = calc(LyrNum(xtalId.getLayer()),
                    xtalId.getColumn(),
                    FaceNum((idents::CalXtalId::XtalFace)(xtalId.getFace())), 
                    diode);
    }
   
    tDiodeIdx(const LyrNum lyr, const ColNum col, const FaceNum face, const DiodeNum diode) :
      TwrWideIndex(calc(lyr,col,face,diode)) {}
     
    tDiodeIdx(const tXtalIdx xtal, const FaceNum face, const DiodeNum diode) :
      TwrWideIndex(calc(xtal.getLyr(),xtal.getCol(),face,diode)) {}
       
    tDiodeIdx(const tXtalIdx xtal, const XtalDiode xDiode) :
      TwrWideIndex(calc(xtal.getLyr(), xtal.getCol(),
                        xDiode.getFace(), xDiode.getDiode())) {}
         
    tDiodeIdx(const tFaceIdx face, const DiodeNum diode) :
      TwrWideIndex(calc(face.getLyr(),face.getCol(),face.getFace(),diode)) {}
           
    tDiodeIdx() : TwrWideIndex() {}
             
    static const unsigned N_VALS = tFaceIdx::N_VALS*DiodeNum::N_VALS;
               
    LyrNum getLyr()   const {return LyrNum((m_data)/LYR_BASE);}
    ColNum getCol()   const {return (m_data%LYR_BASE)/COL_BASE;}
    FaceNum getFace()  const {return FaceNum((idents::CalXtalId::XtalFace)((m_data%COL_BASE)/FACE_BASE));}
    DiodeNum getDiode() const {return DiodeNum(m_data%FACE_BASE);}
                       
                         
    tXtalIdx getTXtalIdx() const {return tXtalIdx(getLyr(),
                                                  getCol());}
                           
    tFaceIdx getTFaceIdx() const {return tFaceIdx(getLyr(),
                                                  getCol(),
                                                  getFace());}

    bool isValid() const {return m_data < N_VALS;}
  private:
    static unsigned calc(const LyrNum lyr, const ColNum col, const FaceNum face, const DiodeNum diode) {
      return lyr.val()*LYR_BASE + col.val()*COL_BASE + face.val()*FACE_BASE + diode.val();
    }
    static const unsigned FACE_BASE = DiodeNum::N_VALS;
    static const unsigned COL_BASE  = FACE_BASE*FaceNum::N_VALS;
    static const unsigned LYR_BASE  = COL_BASE*ColNum::N_VALS;
  };

  /// index all adc channels in single GLAST Cal module
  class tRngIdx : public TwrWideIndex {
  public:
    explicit tRngIdx(const idents::CalXtalId xtalId) {
      if (!xtalId.validFace() || ! xtalId.validRange())
        throw invalid_argument("tRngIdx requires valid face and range info in xtalId.  Programmer error");
      m_data = calc(LyrNum(xtalId.getLayer()),
                    xtalId.getColumn(),
                    FaceNum((idents::CalXtalId::XtalFace)xtalId.getFace()), 
                    RngNum(xtalId.getRange()));
    }
   
    tRngIdx(const LyrNum lyr, const ColNum col, const FaceNum face, const RngNum rng) :
      TwrWideIndex(calc(lyr,col,face,rng)) 
    {}
     
    tRngIdx(const tXtalIdx xtal, const FaceNum face, const RngNum rng) :
      TwrWideIndex(calc(xtal.getLyr(),xtal.getCol(),face,rng)) {}
       
    tRngIdx(const tXtalIdx xtal, const XtalRng xRng) :
      TwrWideIndex(calc(xtal.getLyr(), xtal.getCol(),
                        xRng.getFace(),xRng.getRng())) {}
         
    tRngIdx(const tFaceIdx faceIdx, const RngNum rng) :
      TwrWideIndex(calc(faceIdx.getLyr(),faceIdx.getCol(),faceIdx.getFace(),rng)) {}
           
    tRngIdx() : TwrWideIndex() {}
             
    LyrNum getLyr()  const {return LyrNum((m_data)/LYR_BASE);}
    ColNum getCol()  const {return (m_data%LYR_BASE)/COL_BASE;}
    FaceNum getFace() const {return FaceNum((idents::CalXtalId::XtalFace)((m_data%COL_BASE)/FACE_BASE));}
    RngNum getRng()  const {return RngNum(m_data%FACE_BASE);}
                       
    static const unsigned N_VALS = tFaceIdx::N_VALS*RngNum::N_VALS;
                             
                               
    tXtalIdx getTXtalIdx() const {return tXtalIdx(getLyr(),
                                                  getCol());}
                                 
    tFaceIdx getTFaceIdx() const {return tFaceIdx(getLyr(),
                                                  getCol(),
                                                  getFace());}
                                   
    bool isValid() const {return m_data < N_VALS;}
  private:
    static unsigned calc(const LyrNum lyr, const ColNum col, const FaceNum face, const RngNum rng) {
      return lyr.val()*LYR_BASE + col.val()*COL_BASE + face.val()*FACE_BASE + rng.val();
    }
    static const unsigned FACE_BASE = RngNum::N_VALS;
    static const unsigned COL_BASE  = FACE_BASE*FaceNum::N_VALS;
    static const unsigned LYR_BASE  = COL_BASE*ColNum::N_VALS;
  };

  /** Id all Cal components of same type within full GLAST LAT

      internal integer storage is contiguous and can be used as an array Index
  */

  class LATWideIndex {
  public:  
    /// prefix ++ operator
    LATWideIndex operator++() {
      LATWideIndex tmp(*this);
      m_data++;
      return tmp;
    }

    /// postfix ++ operator
    LATWideIndex& operator++(int) {
      m_data++; 
      return *this;
    }
    
    unsigned val() const {return m_data;}

    bool operator==(const LATWideIndex that) const {
      return m_data == that.m_data;
    }
    
    bool operator<=(const LATWideIndex that) const {
      return m_data <= that.m_data;
    }
    
    bool operator!=(const LATWideIndex that) const {
      return m_data != that.m_data;
    }
    
    bool operator< (const LATWideIndex that) const {
      return m_data <  that.m_data;
    }
    
  protected:
    explicit LATWideIndex(const unsigned val=0) : m_data(val) {}
    unsigned m_data;
  };
  
  /// index class for all Cal layers in GLAST LAT
  class LyrIdx : public LATWideIndex {
  public:
    explicit LyrIdx(const idents::CalXtalId xtalId) :
      LATWideIndex(calc(xtalId.getTower(),
                        LyrNum(xtalId.getLayer()))) {}
    LyrIdx(const TwrNum twr, const LyrNum lyr) :
      LATWideIndex(calc(twr,lyr)) {}
    LyrIdx() : LATWideIndex() {}

    static const unsigned N_VALS  = 
      TwrNum::N_VALS*LyrNum::N_VALS;

    TwrNum getTwr() const {return m_data/TWR_BASE;}
    LyrNum getLyr() const {return LyrNum(m_data%TWR_BASE);}

    bool isValid() const {return m_data < N_VALS;}
    
    std::string toStr() const;

  private:
    static unsigned calc(const TwrNum twr, const LyrNum lyr) {
      return twr.val()*TWR_BASE + lyr.val();
    }
    static const unsigned TWR_BASE  = LyrNum::N_VALS;
  };

  /// index class for all Cal crystals in GLAST LAT
  class XtalIdx : public LATWideIndex {
  public:
    explicit XtalIdx(const idents::CalXtalId xtalId) :
      LATWideIndex(calc(xtalId.getTower(),
                        LyrNum(xtalId.getLayer()),
                        xtalId.getColumn())) {}

    /// inverse of toStr() method
    explicit XtalIdx(const std::string &str);
    
    XtalIdx(const TwrNum twr, const LyrNum lyr, const ColNum col) :
      LATWideIndex(calc(twr,lyr,col)) {}

    XtalIdx(const TwrNum twr,
            const tXtalIdx twrXtalIdx) :
      LATWideIndex(calc(twr,twrXtalIdx.getLyr(), twrXtalIdx.getCol())) {}
                    
    
    XtalIdx() : LATWideIndex() {}

    idents::CalXtalId getCalXtalId() const {
      return idents::CalXtalId(getTwr().val(),
                               getLyr().val(),
                               getCol().val());
    }
    
    tXtalIdx getTXtalIdx() const {
      return tXtalIdx(getLyr(),
                      getCol());
    }

    static const unsigned N_VALS  = 
      TwrNum::N_VALS*LyrNum::N_VALS*ColNum::N_VALS;

    TwrNum getTwr() const {return m_data/TWR_BASE;}
    LyrNum getLyr() const {return LyrNum((m_data%TWR_BASE)/LYR_BASE);}
    ColNum getCol() const {return m_data%LYR_BASE;}

    bool isValid() const {return m_data < N_VALS;}

    std::string toStr() const;

    /// set internal value, make sure you have the right encoding!
    void setVal(const unsigned  val) {
      m_data = val;
    }

  private:
    static unsigned calc(const TwrNum twr, const LyrNum lyr, const ColNum col) {
      return twr.val()*TWR_BASE + lyr.val()*LYR_BASE + col.val();
    }
    static const unsigned LYR_BASE  = ColNum::N_VALS;
    static const unsigned TWR_BASE  = LYR_BASE*LyrNum::N_VALS;
    static const unsigned short N_FIELDS = 3;
  };
  
  /// index class for all Cal xtal faces in GLAST LAT
  class FaceIdx : public LATWideIndex {
  public:
    FaceIdx() : LATWideIndex() {}

    /// construct object from raw index value (make sure you know the
    /// indexing scheme)
    explicit FaceIdx(const unsigned idx) : LATWideIndex(idx) {}

    explicit FaceIdx(const idents::CalXtalId faceId) {
      if (!faceId.validFace())
        throw invalid_argument("FaceIdx requires valid face info in xtalId"
                               ".  Programmer error");

      
      m_data = calc(faceId.getTower(), 
                    LyrNum(faceId.getLayer()),
                    faceId.getColumn(),
                    FaceNum((idents::CalXtalId::XtalFace)faceId.getFace()));
    }

    FaceIdx(const TwrNum twr, const LyrNum lyr, const ColNum col, const FaceNum face) :
      LATWideIndex(calc(twr,lyr,col,face)) {}
    
    FaceIdx(const XtalIdx xtal, const FaceNum face) {
      m_data = xtal.val()*COL_BASE + face.val();    
    }

    idents::CalXtalId getCalXtalId() const {
      return idents::CalXtalId(getTwr().val(),
                               getLyr().val(),
                               getCol().val(),
                               getFace().val());
    }

    static const unsigned N_VALS = XtalIdx::N_VALS*FaceNum::N_VALS;

    XtalIdx getXtalIdx() const {return XtalIdx(getTwr(),
                                               getLyr(),
                                               getCol());}

    TwrNum getTwr()  const {return m_data/TWR_BASE;}
    LyrNum getLyr()  const {return LyrNum((m_data%TWR_BASE)/LYR_BASE);}
    ColNum getCol()  const {return (m_data%LYR_BASE)/COL_BASE;}
    FaceNum getFace() const {return FaceNum((idents::CalXtalId::XtalFace)(m_data%COL_BASE));}

    
    bool isValid() const {return m_data < N_VALS;}

    std::string toStr() const;

  private:
    static unsigned calc(const TwrNum twr, const LyrNum lyr, const ColNum col, const FaceNum face) {
      return twr.val()*TWR_BASE + lyr.val()*LYR_BASE + col.val()*COL_BASE + face.val();
    }
    static const unsigned short COL_BASE  = FaceNum::N_VALS;
    static const unsigned short LYR_BASE  = COL_BASE*ColNum::N_VALS;
    static const unsigned short TWR_BASE  = LYR_BASE*LyrNum::N_VALS;
  };


  /// index class for all Cal photo-diodes in GLAST LAT.
  class DiodeIdx : public LATWideIndex {
  public:
    DiodeIdx(const TwrNum twr, const LyrNum lyr, const ColNum col, const FaceNum face, const DiodeNum diode) :
      LATWideIndex(calc(twr,lyr,col,face,diode)) {}

    DiodeIdx(const XtalIdx xtal, const FaceNum face, const DiodeNum diode) {
      m_data = xtal.val()*COL_BASE + face.val()*FACE_BASE + diode.val();
    }

    DiodeIdx(const XtalIdx xtal, const XtalDiode xDiode) {
      m_data = xtal.val()*COL_BASE + xDiode.val();
    }

    DiodeIdx(const FaceIdx face, const DiodeNum diode) {
      m_data = face.val()*FACE_BASE + diode.val();
    }
    
    DiodeIdx() : LATWideIndex() {}

    /// construct object from raw index value (make sure you know the
    /// indexing scheme)
    explicit DiodeIdx(const unsigned idx) : LATWideIndex(idx) {}

    static const unsigned N_VALS = FaceIdx::N_VALS*DiodeNum::N_VALS;

    TwrNum getTwr()     const {return TwrNum(m_data/TWR_BASE);}
    LyrNum getLyr()     const {return LyrNum((m_data%TWR_BASE)/LYR_BASE);}
    ColNum getCol()     const {return ColNum((m_data%LYR_BASE)/COL_BASE);}
    FaceNum getFace()   const {return FaceNum((idents::CalXtalId::XtalFace)((m_data%COL_BASE)/FACE_BASE));}
    DiodeNum getDiode() const {return DiodeNum(m_data%FACE_BASE);}


    XtalIdx getXtalIdx() const {return XtalIdx(getTwr(),
                                               getLyr(),
                                               getCol());}

    FaceIdx getFaceIdx() const {return FaceIdx(getTwr(),
                                               getLyr(),
                                               getCol(),
                                               getFace());}

    tDiodeIdx getTDiodeIdx() const { return tDiodeIdx(getLyr(),
                                                      getCol(),
                                                      getFace(),
                                                      getDiode());}
    
    bool isValid() const {return m_data < N_VALS;}

    std::string toStr() const;

    /// set internal raw value
    /// \warning (make sure you are using correctly encoded integer)
    void setVal(const unsigned v) {m_data = v;}
  private:
    static unsigned calc(const TwrNum twr, const LyrNum lyr, const ColNum col, 
                         FaceNum face, DiodeNum diode) {
      return twr.val()*TWR_BASE + lyr.val()*LYR_BASE + col.val()*COL_BASE + face.val()*FACE_BASE + diode.val();
    }
    static const unsigned FACE_BASE = DiodeNum::N_VALS;
    static const unsigned COL_BASE  = FACE_BASE*FaceNum::N_VALS;
    static const unsigned LYR_BASE  = COL_BASE*ColNum::N_VALS;
    static const unsigned TWR_BASE  = LYR_BASE*LyrNum::N_VALS;
  };


  /// index class for all Cal ADC channels in GLAST LAT
  class RngIdx : public LATWideIndex {
  public:
    RngIdx(const TwrNum twr, const LyrNum lyr, const ColNum col, const FaceNum face, const RngNum rng) :
      LATWideIndex(calc(twr,lyr,col,face,rng)) 
    {}
    
    RngIdx(const XtalIdx xtal,const  FaceNum face, const RngNum rng) {
      m_data = xtal.val()*COL_BASE + face.val()*FACE_BASE + rng.val();
    }

    RngIdx(const XtalIdx xtal, const XtalRng xRng) {
      m_data = xtal.val()*COL_BASE + xRng.val();
    }

    RngIdx(const FaceIdx faceIdx, const RngNum rng) {
      m_data = faceIdx.val()*FACE_BASE + rng.val();
    }
    
    RngIdx() : LATWideIndex() {}

    /// construct object from raw index value (make sure you know the
    /// indexing scheme)
    explicit RngIdx(const unsigned idx) : LATWideIndex(idx) {}

    idents::CalXtalId getCalXtalId() const {
      return idents::CalXtalId(getTwr().val(),
                               getLyr().val(),
                               getCol().val(),
                               getFace().val(),
                               getRng().val());
    }
    
    TwrNum getTwr()  const {return TwrNum(m_data/TWR_BASE);}
    LyrNum getLyr()  const {return LyrNum((m_data%TWR_BASE)/LYR_BASE);}
    ColNum getCol()  const {return ColNum((m_data%LYR_BASE)/COL_BASE);}
    FaceNum getFace() const {return FaceNum((idents::CalXtalId::XtalFace)((m_data%COL_BASE)/FACE_BASE));}
    RngNum getRng()  const {return RngNum(m_data%FACE_BASE);}

    static const unsigned N_VALS = FaceIdx::N_VALS*RngNum::N_VALS;

    
    XtalIdx getXtalIdx() const {return XtalIdx(getTwr(),
                                               getLyr(),
                                               getCol());}

    FaceIdx getFaceIdx() const {return FaceIdx(getTwr(),
                                               getLyr(),
                                               getCol(),
                                               getFace());}

    DiodeIdx getDiodeIdx() const {return DiodeIdx(getTwr(),
                                                  getLyr(),
                                                  getCol(),
                                                  getFace(),
                                                  getRng().getDiode());}

    bool isValid() const {return m_data < N_VALS;}

    std::string toStr() const;

  private:
    static unsigned calc(const TwrNum twr, const LyrNum lyr, const ColNum col, const FaceNum face, const RngNum rng) {
      return twr.val()*TWR_BASE + lyr.val()*LYR_BASE + col.val()*COL_BASE + face.val()*FACE_BASE + rng.val();
    }
    static const unsigned FACE_BASE = RngNum::N_VALS;
    static const unsigned COL_BASE  = FACE_BASE*FaceNum::N_VALS;
    static const unsigned LYR_BASE  = COL_BASE*ColNum::N_VALS;
    static const unsigned TWR_BASE  = LYR_BASE*LyrNum::N_VALS; 
  }; 

  /** Volume ID field identifiers.  shouldn't be here, but there is no
      global def as far as i can tell.
  */
  enum {fLATObjects, fTowerY, fTowerX, fTowerObjects, fLayer,
        fMeasure, fCALXtal,fCellCmp, fSegment};
   
  /// \note inline b/c func def is outside of class body (due to use of forward declaration)
  inline RngNum DiodeNum::getX8Rng() const {return RngNum(m_data*2);}

  /// \note inline b/c func def is outside of class body (due to use of forward declaration)
  inline RngNum DiodeNum::getX1Rng() const {return RngNum(getX8Rng().val()+1);}


  /// \note inline b/c func def is outside of class body (due to use of forward declaration)
  inline LyrNum::LyrNum(const DirNum dir, const GCRCNum gcrc): 
    SimpleId(gcrc.val()*DirNum::N_VALS + dir.val()) {}

  /// \note inline b/c func def is outside of class body (due to use of forward declaration)
  inline DirNum LyrNum::getDir() const {
    return (m_data%2 == 0) ? X_DIR : Y_DIR;}


  /** \brief id for 4 classes of Cal xtal asymmetry 

      4 asymmetry classes are
      - LL = large diode on both faces
      - LS = large diode on positive face, small diode on neg
      - SL = sm. diode on pos. face, large diode on neg
      - SS = small diode on both faces
      
  */

  class AsymType : public SimpleId {
  public:
    AsymType(const DiodeNum posDiode, const DiodeNum negDiode) :
      SimpleId(posDiode.val()*2 + negDiode.val()) {};

    /// instantiate class from raw interan index (better have correct index!)
    explicit AsymType(const unsigned short val=0) :
      SimpleId(val) 
    {}
    
    DiodeNum getDiode(const FaceNum face) const {
      switch ((idents::CalXtalId::XtalFace)face.val()) {
      case idents::CalXtalId::POS:
        return DiodeNum(m_data/2);
      case idents::CalXtalId::NEG:
        return DiodeNum(m_data%2);
      default:
        throw invalid_argument("Bad FaceNum value. Programmer error");
      }
    }

    static const unsigned short N_VALS = 4; 

    bool isValid() const {return m_data < N_VALS;}

    const std::string &toStr() const;


    bool operator==(const AsymType that) const {return m_data == that.m_data;}
    bool operator!=(const AsymType that) const {return m_data != that.m_data;}

  };
  const AsymType ASYM_LL(LRG_DIODE, LRG_DIODE);
  const AsymType ASYM_LS(LRG_DIODE, SM_DIODE);
  const AsymType ASYM_SL(SM_DIODE, LRG_DIODE);
  const AsymType ASYM_SS(SM_DIODE, SM_DIODE);

  /// index class for all Cal GCFE's
  class GCFEIdx : public LATWideIndex {
  public:
    /// \param val should be the result of GCFEIdx::val() method to ensure proper encoding.
    GCFEIdx(const unsigned val=0) : LATWideIndex(val) {}

    /// inverse of toStr() method
    explicit GCFEIdx(const std::string &str);
    
    GCFEIdx(const GTEMNum tem,
            const GCCCNum ccc,
            const GCRCNum crc,
            const GCFENum cfe) :
      LATWideIndex(calc(tem, ccc, crc, cfe)) {}

    static const unsigned N_VALS  = 
      GTEMNum::N_VALS*GCCCNum::N_VALS*GCRCNum::N_VALS*GCFENum::N_VALS;

    GTEMNum getGTEM() const {return m_data/TEM_BASE;}
    GCCCNum getGCCC() const {return GCCCNum((m_data%TEM_BASE)/CCC_BASE);}
    GCRCNum getGCRC() const {return GCRCNum((m_data%CCC_BASE)/CRC_BASE);}
    GCFENum getGCFE() const {return m_data%CRC_BASE;}

    bool isValid() const {return m_data < N_VALS;}

    std::string toStr() const;

  private:
    static unsigned calc(const GTEMNum tem,
                         const GCCCNum ccc,
                         const GCRCNum crc,
                         const GCFENum cfe) {
      return tem.val()*TEM_BASE + ccc.val()*CCC_BASE + crc.val()*CRC_BASE + cfe.val()*CFE_BASE;
    }

    static const unsigned CFE_BASE = 1;
    static const unsigned CRC_BASE = GCFENum::N_VALS;
    static const unsigned CCC_BASE = CRC_BASE*GCRCNum::N_VALS;
    static const unsigned TEM_BASE = CCC_BASE*GCCCNum::N_VALS;
    static const unsigned short N_FIELDS = 4;

  };

  /// index class for all Cal GCRC's
  class GCRCIdx : public LATWideIndex {
  public:
    /// \param val should be the result of GCRCIdx::val() method to ensure proper encoding.
    /// \note explicit as row indexing is ambiguous
    explicit GCRCIdx(const unsigned val=0) : LATWideIndex(val) {}

    /// inverse of toStr() method
    explicit GCRCIdx(const std::string &str);
    
    GCRCIdx(const GTEMNum tem,
            const GCCCNum ccc,
            const GCRCNum crc) :
      LATWideIndex(calc(tem, ccc, crc)) {}

    static const unsigned N_VALS  = 
      GTEMNum::N_VALS*GCCCNum::N_VALS*GCRCNum::N_VALS;

    GTEMNum getGTEM() const {return m_data/TEM_BASE;}
    GCCCNum getGCCC() const {return GCCCNum((m_data%TEM_BASE)/CCC_BASE);}
    GCRCNum getGCRC() const {return GCRCNum(m_data%CCC_BASE);}

    bool isValid() const {return m_data < N_VALS;}

    std::string toStr() const;

  private:
    static unsigned calc(const GTEMNum tem,
                         const GCCCNum ccc,
                         const GCRCNum crc) {
      return tem.val()*TEM_BASE + ccc.val()*CCC_BASE + crc.val()*CRC_BASE;
    }

    static const unsigned CRC_BASE = 1;
    static const unsigned CCC_BASE = CRC_BASE*GCRCNum::N_VALS;
    static const unsigned TEM_BASE = CCC_BASE*GCCCNum::N_VALS;
    static const unsigned short N_FIELDS = 3;
  };

  /// index class for all Cal GCCC's
  class GCCCIdx : public LATWideIndex {
  public:
    /// \param val should be the result of GCCCIdx::val() method to ensure proper encoding.
    GCCCIdx(const unsigned val=0) : LATWideIndex(val) {}

    /// inverse of toStr() method
    explicit GCCCIdx(const std::string &str);
    
    GCCCIdx(const GTEMNum tem,
            const GCCCNum ccc) :
      LATWideIndex(calc(tem, ccc)) {}

    static const unsigned N_VALS  = 
      GTEMNum::N_VALS*GCCCNum::N_VALS;

    GTEMNum getGTEM() const {return m_data/TEM_BASE;}
    GCCCNum getGCCC() const {return GCCCNum(m_data%TEM_BASE);}

    bool isValid() const {return m_data < N_VALS;}

    std::string toStr() const;

  private:
    static unsigned calc(const GTEMNum tem,
                         const GCCCNum ccc) {
      return tem.val()*TEM_BASE + ccc.val()*CCC_BASE;
    }

    static const unsigned CCC_BASE = 1;
    static const unsigned TEM_BASE = CCC_BASE*GCCCNum::N_VALS;
    static const unsigned short N_FIELDS = 2;
  };
  
  /// Saturation of CAL ADC.
  /// We deliberately choose a value slightly lower than 2**12 - 1 = 4096, 
  static const float CAL_ADC_SATURATION = 4060;
  

}; // namespace CalUtil


#endif // CalDefs_H
