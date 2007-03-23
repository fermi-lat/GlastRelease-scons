#ifndef CalDefs_H
#define CalDefs_H

// LOCAL

// GLAST
#include "idents/CalXtalId.h"

// EXTLIB

// STD
#include <string>
#include <vector>
#include <ostream>
#include <stdexcept>

using namespace std;

/** @file CalDefs.h
    \author Zachary Fewtrell

    \brief enums & array index types for use throughout GLAST Cal software.
*/
namespace CalUtil {
  // forward declarations (sometimes needed)
  class DirNum;
  class RngNum;

  /** \brief Atomic (simple) Cal Geometry Id's

  Id single Cal component w/ in it's immediate container
  */
  class SimpleId {
  public:
    operator unsigned short() const {return m_data;}

    /// prefix ++ operator
    SimpleId operator++() {
      SimpleId tmp(*this);
      m_data++;
      return tmp;
    }

    /// postfix ++ operator
    SimpleId& operator++(int) {
      m_data++; 
      return *this;
    }

    unsigned short val() const {return m_data;}
  protected:
    SimpleId(unsigned short val) : m_data(val) {}
    SimpleId() : m_data(0) {}
    unsigned short m_data;
  };

  /// index class for GLAST Cal module tower bay
  class TwrNum : public SimpleId {
  public:
    TwrNum() : SimpleId() {}
    TwrNum(unsigned short val) : SimpleId(val) {}
    TwrNum(unsigned short tRow, unsigned short tCol) : 
      SimpleId(tRow*N_COLS + tCol) {}

    unsigned short  getRow() const {return m_data/N_COLS;}
    unsigned short  getCol() const {return m_data%N_COLS;}

    static const unsigned short N_VALS = 16; 
    bool isValid() const {return m_data < N_VALS;}
            
    static const unsigned short N_COLS = 4;
    static const unsigned short N_ROWS = 4;

    bool operator==(const TwrNum &that) const {return m_data == that.m_data;}
    bool operator!=(const TwrNum &that) const {return m_data != that.m_data;}
    bool operator>=(const TwrNum &that) const {return m_data >= that.m_data;}
    bool operator<=(const TwrNum &that) const {return m_data <= that.m_data;}
    bool operator>(const TwrNum &that) const {return m_data > that.m_data;}
    bool operator<(const TwrNum &that) const {return m_data < that.m_data;}
  };

  /// index class for GLAST Cal GCRC (layer 0-3 on a single cal front end board).
  class GCRCNum : public SimpleId {
  public:
    GCRCNum() : SimpleId() {}
    GCRCNum(unsigned short val) : SimpleId(val) {}

    static const unsigned short N_VALS=4;
    bool isValid() const {return m_data < N_VALS;}

    bool operator==(const GCRCNum &that) const {return m_data == that.m_data;}
    bool operator!=(const GCRCNum &that) const {return m_data != that.m_data;}
    bool operator>=(const GCRCNum &that) const {return m_data >= that.m_data;}
    bool operator<=(const GCRCNum &that) const {return m_data <= that.m_data;}
    bool operator>(const GCRCNum &that) const {return m_data > that.m_data;}
    bool operator<(const GCRCNum &that) const {return m_data < that.m_data;}
  };

  /// index class for GLAST Cal xtal layer w/in Cal module
  class LyrNum : public SimpleId {
  public:
    LyrNum() : SimpleId() {}
    LyrNum(unsigned short val) : SimpleId(val) {}
    inline LyrNum(DirNum dir, GCRCNum gcrc);
    
    static const unsigned short N_VALS=8;
    bool isValid() const {return m_data < N_VALS;}

    // 0 is 'X' Direction
    DirNum getDir() const;
    GCRCNum getGCRC() const {return GCRCNum(m_data/2);}

    bool operator==(const LyrNum &that) const {return m_data == that.m_data;}
    bool operator!=(const LyrNum &that) const {return m_data != that.m_data;}
    bool operator>=(const LyrNum &that) const {return m_data >= that.m_data;}
    bool operator<=(const LyrNum &that) const {return m_data <= that.m_data;}
    bool operator>(const LyrNum &that) const {return m_data > that.m_data;}
    bool operator<(const LyrNum &that) const {return m_data < that.m_data;}

  };

  /// index class for GLAST Cal xtal direction ('X' or 'Y')
  class DirNum : public SimpleId {
  public:
    DirNum() : SimpleId() {}
    DirNum(unsigned short val) : SimpleId(val) {}

    static const unsigned short N_VALS=2;
    bool isValid() const {return m_data < N_VALS;}

    bool operator==(const DirNum &that) const {return m_data == that.m_data;}
    bool operator!=(const DirNum &that) const {return m_data != that.m_data;}

	const std::string &toStr() const;

  };
  const DirNum X_DIR(0);
  const DirNum Y_DIR(1);

  /// index class for GLAST Cal xtal column w/in Cal xtal layer
  class ColNum : public SimpleId {
  public:
    ColNum() : SimpleId() {}
    ColNum(unsigned short val) : SimpleId(val) {}

    static const unsigned short N_VALS=12;
    bool isValid() const {return m_data < N_VALS;}

    bool operator==(const ColNum &that) const {return m_data == that.m_data;}
    bool operator!=(const ColNum &that) const {return m_data != that.m_data;}
    bool operator>=(const ColNum &that) const {return m_data >= that.m_data;}
    bool operator<=(const ColNum &that) const {return m_data <= that.m_data;}
    bool operator>(const ColNum &that) const {return m_data > that.m_data;}
    bool operator<(const ColNum &that) const {return m_data < that.m_data;}

  };

  /// index class for GLAST Cal xtal face ('POS' or 'NEG')
  class FaceNum : public SimpleId {
  public:
    FaceNum() : SimpleId() {}
    FaceNum(unsigned short val) : SimpleId(val) {}

    static const vector<string> MNEM;
    static const unsigned short N_VALS=2;    
    bool isValid() const {return m_data < N_VALS;}
	const std::string &toStr() const;
    
    /// allow quick conversion to idents::CalXtalId::XtalFace since internal
    /// storage is the same
    operator idents::CalXtalId::XtalFace() {
      return (idents::CalXtalId::XtalFace)m_data;
    }

    bool operator==(const FaceNum &that) const {return m_data == that.m_data;}
    bool operator!=(const FaceNum &that) const {return m_data != that.m_data;}

  };
  const FaceNum POS_FACE(idents::CalXtalId::POS);
  const FaceNum NEG_FACE(idents::CalXtalId::NEG);

  /// index class for GLAST Cal photo-diode on single xtal face 
  /// ('LARGE' or 'SMALL')
  class DiodeNum : public SimpleId {
  public:
    DiodeNum(unsigned short val) : SimpleId(val) {}
    DiodeNum() : SimpleId() {}
    
    inline RngNum getX8Rng() const;
    inline RngNum getX1Rng() const;

    static const vector<string> MNEM;
    static const unsigned short N_VALS = 2; 
    bool isValid() const {return m_data < N_VALS;}

	const std::string &toStr() const;

    bool operator==(const DiodeNum &that) const {return m_data == that.m_data;}
    bool operator!=(const DiodeNum &that) const {return m_data != that.m_data;}

  };
  const DiodeNum LRG_DIODE(idents::CalXtalId::LARGE);
  const DiodeNum SM_DIODE (idents::CalXtalId::SMALL);

  /// THX (Track & Hold Multiplier) can be either X8 or X1
  class THXNum : public SimpleId {
  public:
    THXNum(unsigned short val) : SimpleId(val) {}
    THXNum() : SimpleId() {}

    static const vector<string> MNEM;

    bool isValid() const {return m_data < N_VALS;}
    static const unsigned short N_VALS=2;
	const std::string &toStr() const;

    bool operator==(const THXNum &that) const {return m_data == that.m_data;}
    bool operator!=(const THXNum &that) const {return m_data != that.m_data;}

  };
  const THXNum THX8(0);
  const THXNum THX1(1);

  /// index class for GLAST Cal ADC range (LEX8 -> HEX1)
  class RngNum : public SimpleId {
  public:
    RngNum(unsigned short val) : SimpleId(val) {}
    RngNum(DiodeNum diode, THXNum thx) :
      SimpleId(diode.val()*2 + thx) {}
    RngNum() : SimpleId() {}

    static const vector<string> MNEM;
    static const unsigned short N_VALS=4;
	const std::string &toStr() const;
    
    DiodeNum getDiode() const {
      using idents::CalXtalId;
      return idents::CalXtalId::rangeToDiode((idents::CalXtalId::AdcRange)m_data);
    }
    bool isValid() const {return m_data < N_VALS;}

    operator idents::CalXtalId::AdcRange() {
      return (idents::CalXtalId::AdcRange)m_data;
    }

	THXNum getTHX() {return m_data%2;}

    bool operator==(const RngNum &that) const {return m_data == that.m_data;}
    bool operator!=(const RngNum &that) const {return m_data != that.m_data;}
    bool operator>=(const RngNum &that) const {return m_data >= that.m_data;}
    bool operator<=(const RngNum &that) const {return m_data <= that.m_data;}
    bool operator>(const RngNum &that) const {return m_data > that.m_data;}
    bool operator<(const RngNum &that) const {return m_data < that.m_data;}

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
    XtalWideIndex& operator++(int) {
      m_data++; 
      return *this;
    }

    bool operator==(const XtalWideIndex &that) const {
      return m_data == that.m_data;}
    bool operator<=(const XtalWideIndex &that) const {
      return m_data <= that.m_data;}
    bool operator!=(const XtalWideIndex &that) const {
      return m_data != that.m_data;}
    bool operator< (const XtalWideIndex &that) const {
      return m_data <  that.m_data;}

    unsigned val() const {return m_data;}

  protected:
    XtalWideIndex(unsigned short val) : m_data(val) {}
    XtalWideIndex() : m_data(0) {}
    unsigned short m_data;
  };

  /// index class for all 4 photo diodes in one Cal xtal.
  class XtalDiode : public XtalWideIndex {
  public: 
    XtalDiode(FaceNum face, DiodeNum diode) :
      XtalWideIndex(face.val()*FACE_BASE + diode.val()) {}
    XtalDiode() : XtalWideIndex() {}

    DiodeNum getDiode() const {return m_data%FACE_BASE;}
    FaceNum getFace()  const {return m_data/FACE_BASE;}

    static const unsigned short N_VALS = FaceNum::N_VALS*DiodeNum::N_VALS;
    bool isValid() const {return m_data < N_VALS;}
  protected:
    static const unsigned short FACE_BASE = DiodeNum::N_VALS;
  };

  /// index class for all 8 ADC readouts on one Cal xtal.
  class XtalRng : public XtalWideIndex {
  public:
    XtalRng(FaceNum face, RngNum rng) :
      XtalWideIndex(face.val()*FACE_BASE + rng.val()) {};
    XtalRng() : XtalWideIndex() {}

    FaceNum getFace() const {return m_data/FACE_BASE;}
    RngNum  getRng() const {return m_data%FACE_BASE;}

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
       
    bool operator==(const TwrWideIndex &that) const {return m_data == that.m_data;}
    bool operator<=(const TwrWideIndex &that) const {return m_data <= that.m_data;}
    bool operator!=(const TwrWideIndex &that) const {return m_data != that.m_data;}
    bool operator< (const TwrWideIndex &that) const {return m_data <  that.m_data;}
    //TwrWideIndex& operator= (const TwrWideIndex &that) {m_data = that.m_data;}
  protected:
    TwrWideIndex(unsigned val) : m_data(val) {}
    TwrWideIndex() : m_data(0) {}
    unsigned m_data;
  };

  /// index all crytals in single cal module
  class tXtalIdx : public TwrWideIndex {
  public:
    tXtalIdx(const idents::CalXtalId &xtal) :
      TwrWideIndex(calc(xtal.getLayer(),
                        xtal.getColumn())) {}
   
    tXtalIdx(LyrNum lyr, ColNum col) :
      TwrWideIndex(calc(lyr, col)) {}
       
    tXtalIdx() : TwrWideIndex() {}
       
    idents::CalXtalId getCalXtalId() const {
      return idents::CalXtalId(0,//twr
                               getLyr(),
                               getCol());
    }
    static const unsigned N_VALS  = LyrNum::N_VALS*ColNum::N_VALS;
           
    LyrNum getLyr() const {return (m_data)/LYR_BASE;}
    ColNum getCol() const {return m_data%LYR_BASE;}
               
    bool isValid() const {return m_data < N_VALS;}
                   
  private:
    static unsigned calc(LyrNum lyr, ColNum col) {
      return lyr*LYR_BASE + col;
    }
    static const unsigned LYR_BASE  = ColNum::N_VALS;
  };

  /// index all crystal faces in single cal module
  class tFaceIdx : public TwrWideIndex {
  public:
    tFaceIdx() : TwrWideIndex() {}
   
    tFaceIdx(const idents::CalXtalId &xtalId) {
      if (!xtalId.validFace())
        throw invalid_argument("tFaceIdx requires valid face info in xtalId."
                               "  Programmer error");
                           
      m_data = calc(xtalId.getLayer(),
                    xtalId.getColumn(),
                    xtalId.getFace());
    }
     
    tFaceIdx(LyrNum lyr, ColNum col, FaceNum face) :
      TwrWideIndex(calc(lyr,col,face)) {}
       
    tFaceIdx(tXtalIdx xtal, FaceNum face) :
      TwrWideIndex(calc(xtal.getLyr(), xtal.getCol(), face)) {}
         
    idents::CalXtalId getCalXtalId() const {
      return idents::CalXtalId(0,//twr
                               getLyr(),
                               getCol(),
                               getFace().val());
    }
           
    static const unsigned N_VALS = tXtalIdx::N_VALS*FaceNum::N_VALS;
             
    tXtalIdx getTXtalIdx() const {return tXtalIdx(getLyr(),
                                                  getCol());}
               
    LyrNum getLyr()  const {return (m_data)/LYR_BASE;}
    ColNum getCol()  const {return (m_data%LYR_BASE)/COL_BASE;}
    FaceNum getFace() const {return m_data%COL_BASE;}
                     
    void setFace(FaceNum face) {
      m_data = calc(getLyr(),getCol(),face);
    }
                       
                         
    bool isValid() const {return m_data < N_VALS;}
  private:
    static unsigned calc(LyrNum lyr, ColNum col, FaceNum face) {
      return lyr.val()*LYR_BASE + col.val()*COL_BASE + face.val();
    }
    static const unsigned short COL_BASE  = FaceNum::N_VALS;
    static const unsigned short LYR_BASE  = COL_BASE*ColNum::N_VALS;
  };

  /// index all crystal diodes in single cal module
  class tDiodeIdx : public TwrWideIndex {
  public:
    tDiodeIdx(const idents::CalXtalId &xtalId, DiodeNum diode) {
      if (!xtalId.validFace())
        throw invalid_argument("tDiodeIdx requires valid face info in xtalId.  Programmer error");
                           
      m_data = calc(xtalId.getLayer(),
                    xtalId.getColumn(),
                    xtalId.getFace(), 
                    diode);
    }
   
    tDiodeIdx(LyrNum lyr, ColNum col, FaceNum face, DiodeNum diode) :
      TwrWideIndex(calc(lyr,col,face,diode)) {}
     
    tDiodeIdx(tXtalIdx xtal, FaceNum face, DiodeNum diode) :
      TwrWideIndex(calc(xtal.getLyr(),xtal.getCol(),face,diode)) {}
       
    tDiodeIdx(tXtalIdx xtal, XtalDiode xDiode) :
      TwrWideIndex(calc(xtal.getLyr(), xtal.getCol(),
                        xDiode.getFace(), xDiode.getDiode())) {}
         
    tDiodeIdx(tFaceIdx face, DiodeNum diode) :
      TwrWideIndex(calc(face.getLyr(),face.getCol(),face.getFace(),diode)) {}
           
    tDiodeIdx() : TwrWideIndex() {}
             
    static const unsigned N_VALS = tFaceIdx::N_VALS*DiodeNum::N_VALS;
               
    LyrNum getLyr()   const {return (m_data)/LYR_BASE;}
    ColNum getCol()   const {return (m_data%LYR_BASE)/COL_BASE;}
    FaceNum getFace()  const {return (m_data%COL_BASE)/FACE_BASE;}
    DiodeNum getDiode() const {return m_data%FACE_BASE;}
                       
                         
    tXtalIdx getTXtalIdx() const {return tXtalIdx(getLyr(),
                                                  getCol());}
                           
    tFaceIdx getTFaceIdx() const {return tFaceIdx(getLyr(),
                                                  getCol(),
                                                  getFace());}

    tDiodeIdx getTDiodeIdx() const {return tDiodeIdx(getLyr(),
                                                    getCol(),
                                                    getFace(),
                                                    getDiode());}
                             
                             
    bool isValid() const {return m_data < N_VALS;}
  private:
    static unsigned calc(LyrNum lyr, ColNum col, FaceNum face, DiodeNum diode) {
      return lyr.val()*LYR_BASE + col.val()*COL_BASE + face.val()*FACE_BASE + diode.val();
    }
    static const unsigned FACE_BASE = DiodeNum::N_VALS;
    static const unsigned COL_BASE  = FACE_BASE*FaceNum::N_VALS;
    static const unsigned LYR_BASE  = COL_BASE*ColNum::N_VALS;
  };

  /// index all adc channels in single GLAST Cal module
  class tRngIdx : public TwrWideIndex {
  public:
    tRngIdx(const idents::CalXtalId &xtalId) {
      if (!xtalId.validFace() || ! xtalId.validRange())
        throw invalid_argument("tRngIdx requires valid face and range info in xtalId.  Programmer error");
      m_data = calc(xtalId.getLayer(),xtalId.getColumn(),xtalId.getFace(), xtalId.getRange());
    }
   
    tRngIdx(LyrNum lyr, ColNum col, FaceNum face, RngNum rng) :
      TwrWideIndex(calc(lyr,col,face,rng)) 
      {}
     
    tRngIdx(tXtalIdx xtal, FaceNum face, RngNum rng) :
      TwrWideIndex(calc(xtal.getLyr(),xtal.getCol(),face,rng)) {}
       
    tRngIdx(tXtalIdx xtal, XtalRng xRng) :
      TwrWideIndex(calc(xtal.getLyr(), xtal.getCol(),
                        xRng.getFace(),xRng.getRng())) {}
         
    tRngIdx(tFaceIdx faceIdx, RngNum rng) :
      TwrWideIndex(calc(faceIdx.getLyr(),faceIdx.getCol(),faceIdx.getFace(),rng)) {}
           
    tRngIdx() : TwrWideIndex() {}
             
    idents::CalXtalId getCalXtalId() const {
      return idents::CalXtalId(getLyr(),
                               getCol(),
                               getFace().val(),
                               getRng().val());
    }
               
    LyrNum getLyr()  const {return (m_data)/LYR_BASE;}
    ColNum getCol()  const {return (m_data%LYR_BASE)/COL_BASE;}
    FaceNum getFace() const {return (m_data%COL_BASE)/FACE_BASE;}
    RngNum getRng()  const {return m_data%FACE_BASE;}
                       
    void setFace(FaceNum face) {
      m_data = calc(getLyr(),getCol(),face,getRng());
    }
    void setRng(RngNum rng) {
      m_data = calc(getLyr(),getCol(),getFace(),rng);
    }
                           
    static const unsigned N_VALS = tFaceIdx::N_VALS*RngNum::N_VALS;
                             
                               
    tXtalIdx getTXtalIdx() const {return tXtalIdx(getLyr(),
                                                  getCol());}
                                 
    tFaceIdx getTFaceIdx() const {return tFaceIdx(getLyr(),
                                                  getCol(),
                                                  getFace());}
                                   
    bool isValid() const {return m_data < N_VALS;}
  private:
    static unsigned calc(LyrNum lyr, ColNum col, FaceNum face, RngNum rng) {
      return lyr*LYR_BASE + col*COL_BASE + face.val()*FACE_BASE + rng.val();
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

	/// direct set of internal integer value 
	void setVal(unsigned val) {m_data = val;}
    
    bool operator==(const LATWideIndex &that) const {
      return m_data == that.m_data;
    }
    bool operator<=(const LATWideIndex &that) const {
      return m_data <= that.m_data;
    }
    bool operator!=(const LATWideIndex &that) const {
      return m_data != that.m_data;
    }
    bool operator< (const LATWideIndex &that) const {
      return m_data <  that.m_data;
    }
    //LATWideIndex& operator= (const LATWideIndex &that) {
    //m_data = that.m_data;}
  protected:
    LATWideIndex(unsigned val) : m_data(val) {}
    LATWideIndex() : m_data(0) {}
    unsigned m_data;
  };
  
  /// index class for all Cal crystals in GLAST LAT
  class LyrIdx : public LATWideIndex {
  public:
    LyrIdx(const idents::CalXtalId &xtalId) :
      LATWideIndex(calc(xtalId.getTower(),
                        xtalId.getLayer())) {}
    LyrIdx(TwrNum twr, LyrNum lyr) :
      LATWideIndex(calc(twr,lyr)) {}
    LyrIdx() : LATWideIndex() {}

    static const unsigned N_VALS  = 
      TwrNum::N_VALS*LyrNum::N_VALS;

    TwrNum getTwr() const {return m_data/TWR_BASE;}
    LyrNum getLyr() const {return m_data%TWR_BASE;}

    bool isValid() const {return m_data < N_VALS;}
    
  private:
    static unsigned calc(TwrNum twr, LyrNum lyr) {
      return twr*TWR_BASE + lyr;
    }
    static const unsigned TWR_BASE  = LyrNum::N_VALS;
  };

  /// index class for all Cal crystals in GLAST LAT
  class XtalIdx : public LATWideIndex {
  public:
    XtalIdx(const idents::CalXtalId &xtalId) :
      LATWideIndex(calc(xtalId.getTower(),
                        xtalId.getLayer(),
                        xtalId.getColumn())) {}
    XtalIdx(TwrNum twr, LyrNum lyr, ColNum col) :
      LATWideIndex(calc(twr,lyr,col)) {}
    XtalIdx() : LATWideIndex() {}

    idents::CalXtalId getCalXtalId() const {
      return idents::CalXtalId(getTwr(),
                               getLyr(),
                               getCol());
    }
    
    tXtalIdx getTXtalIdx() const {
      return tXtalIdx(getLyr(),
                      getCol());
    }

    static const unsigned N_VALS  = 
      TwrNum::N_VALS*LyrNum::N_VALS*ColNum::N_VALS;

    TwrNum getTwr() const {return m_data/TWR_BASE;}
    LyrNum getLyr() const {return (m_data%TWR_BASE)/LYR_BASE;}
    ColNum getCol() const {return m_data%LYR_BASE;}

    bool isValid() const {return m_data < N_VALS;}

	std::string toStr() const;
    
  private:
    static unsigned calc(TwrNum twr, LyrNum lyr, ColNum col) {
      return twr*TWR_BASE + lyr*LYR_BASE + col;
    }
    static const unsigned LYR_BASE  = ColNum::N_VALS;
    static const unsigned TWR_BASE  = LYR_BASE*LyrNum::N_VALS;
  };
  
  /// index class for all Cal xtal faces in GLAST LAT
  class FaceIdx : public LATWideIndex {
  public:
    FaceIdx() : LATWideIndex() {}
    FaceIdx(const idents::CalXtalId &faceId) {
      if (!faceId.validFace())
        throw invalid_argument("FaceIdx requires valid face info in xtalId"
                               ".  Programmer error");

      
      m_data = calc(faceId.getTower(), 
                    faceId.getLayer(),
                    faceId.getColumn(),
                    faceId.getFace());
    }
    FaceIdx(TwrNum twr, LyrNum lyr, ColNum col, FaceNum face) :
      LATWideIndex(calc(twr,lyr,col,face)) {}
    FaceIdx(XtalIdx xtal, FaceNum face) {
      m_data = xtal.val()*COL_BASE + face.val();    
    }

    idents::CalXtalId getCalXtalId() const {
      return idents::CalXtalId(getTwr(),
                               getLyr(),
                               getCol(),
                               getFace().val());
    }

    static const unsigned N_VALS = XtalIdx::N_VALS*FaceNum::N_VALS;

    XtalIdx getXtalIdx() const {return XtalIdx(getTwr(),
                                               getLyr(),
                                               getCol());}

    TwrNum getTwr()  const {return m_data/TWR_BASE;}
    LyrNum getLyr()  const {return (m_data%TWR_BASE)/LYR_BASE;}
    ColNum getCol()  const {return (m_data%LYR_BASE)/COL_BASE;}
    FaceNum getFace() const {return m_data%COL_BASE;}

    
    bool isValid() const {return m_data < N_VALS;}

	std::string toStr() const;
  private:
    static unsigned calc(TwrNum twr, LyrNum lyr, ColNum col, FaceNum face) {
      return twr*TWR_BASE + lyr*LYR_BASE + col*COL_BASE + face.val();
    }
    static const unsigned short COL_BASE  = FaceNum::N_VALS;
    static const unsigned short LYR_BASE  = COL_BASE*ColNum::N_VALS;
    static const unsigned short TWR_BASE  = LYR_BASE*LyrNum::N_VALS;
  };


  /// index class for all Cal photo-diodes in GLAST LAT.
  class DiodeIdx : public LATWideIndex {
  public:
    DiodeIdx(TwrNum twr, LyrNum lyr, ColNum col, FaceNum face, DiodeNum diode) :
      LATWideIndex(calc(twr,lyr,col,face,diode)) {}

    DiodeIdx(XtalIdx xtal, FaceNum face, DiodeNum diode) {
      m_data = xtal.val()*COL_BASE + face.val()*FACE_BASE + diode.val();
    }

    DiodeIdx(XtalIdx xtal, XtalDiode xDiode) {
      m_data = xtal.val()*COL_BASE + xDiode.val();
    }

    DiodeIdx(FaceIdx face, DiodeNum diode) {
      m_data = face.val()*FACE_BASE + diode.val();
    }
    
    DiodeIdx() : LATWideIndex() {}

    static const unsigned N_VALS = FaceIdx::N_VALS*DiodeNum::N_VALS;

    TwrNum getTwr()     const {return m_data/TWR_BASE;}
    LyrNum getLyr()     const {return (m_data%TWR_BASE)/LYR_BASE;}
    ColNum getCol()     const {return (m_data%LYR_BASE)/COL_BASE;}
    FaceNum getFace()   const {return (m_data%COL_BASE)/FACE_BASE;}
    DiodeNum getDiode() const {return m_data%FACE_BASE;}


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
  private:
    static unsigned calc(TwrNum twr, LyrNum lyr, ColNum col, 
                         FaceNum face, DiodeNum diode) {
      return twr*TWR_BASE + lyr*LYR_BASE + col*COL_BASE + face.val()*FACE_BASE + diode.val();
    }
    static const unsigned FACE_BASE = DiodeNum::N_VALS;
    static const unsigned COL_BASE  = FACE_BASE*FaceNum::N_VALS;
    static const unsigned LYR_BASE  = COL_BASE*ColNum::N_VALS;
    static const unsigned TWR_BASE  = LYR_BASE*LyrNum::N_VALS;
  };


  /// index class for all Cal ADC channels in GLAST LAT
  class RngIdx : public LATWideIndex {
  public:
    RngIdx(TwrNum twr, LyrNum lyr, ColNum col, FaceNum face, RngNum rng) :
      LATWideIndex(calc(twr,lyr,col,face,rng)) 
      {}
    
    RngIdx(XtalIdx xtal, FaceNum face, RngNum rng) {
      m_data = xtal.val()*COL_BASE + face.val()*FACE_BASE + rng.val();
    }

    RngIdx(XtalIdx xtal, XtalRng xRng) {
      m_data = xtal.val()*COL_BASE + xRng.val();
    }

    RngIdx(FaceIdx faceIdx, RngNum rng) {
      m_data = faceIdx.val()*FACE_BASE + rng.val();
    }
    
    RngIdx() : LATWideIndex() {}

    idents::CalXtalId getCalXtalId() const {
      return idents::CalXtalId(getTwr(),
                               getLyr(),
                               getCol(),
                               getFace().val(),
                               getRng().val());
    }
    
    TwrNum getTwr()  const {return m_data/TWR_BASE;}
    LyrNum getLyr()  const {return (m_data%TWR_BASE)/LYR_BASE;}
    ColNum getCol()  const {return (m_data%LYR_BASE)/COL_BASE;}
    FaceNum getFace() const {return (m_data%COL_BASE)/FACE_BASE;}
    RngNum getRng()  const {return m_data%FACE_BASE;}

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
    static unsigned calc(TwrNum twr, LyrNum lyr, ColNum col, FaceNum face, RngNum rng) {
      return twr*TWR_BASE + lyr*LYR_BASE + col*COL_BASE + face.val()*FACE_BASE + rng.val();
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
   
  // some functions must be placed after definition of types
  inline RngNum DiodeNum::getX8Rng() const {return m_data*2;}

  inline RngNum DiodeNum::getX1Rng() const {return getX8Rng().val()+1;}

  inline LyrNum::LyrNum(DirNum dir, GCRCNum gcrc): 
    SimpleId((unsigned short)gcrc*DirNum::N_VALS + (unsigned short)dir) {}


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
    AsymType() : SimpleId() {}
    AsymType(DiodeNum posDiode, DiodeNum negDiode) :
      SimpleId(posDiode.val()*2 + negDiode.val()) {};
    
    inline DiodeNum getDiode(FaceNum face) {
      switch ((idents::CalXtalId::XtalFace)face) {
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


    bool operator==(const AsymType &that) const {return m_data == that.m_data;}
    bool operator!=(const AsymType &that) const {return m_data != that.m_data;}

  };
  const AsymType ASYM_LL(LRG_DIODE, LRG_DIODE);
  const AsymType ASYM_LS(LRG_DIODE, SM_DIODE);
  const AsymType ASYM_SL(SM_DIODE, LRG_DIODE);
  const AsymType ASYM_SS(SM_DIODE, SM_DIODE);

}; // namespace CalUtil


#endif // CalDefs_H

