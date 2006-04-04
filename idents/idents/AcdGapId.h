#ifndef idents_AcdGapId_H
#define idents_AcdGapId_H

/** @class AcdGapId
 @brief encapsulate the id for an ACD detector. 

*/

namespace idents {

  class AcdGapId {
    
  private:
    
    enum { ColShift = 0,
	   RowShift = 3,
	   FaceShift = 6,
	   GapShift = 9,
	   TypeShift = 12,
	   ThreeBitMask = 0x7,
	   NibbleMask = 0xF};

  public:
    
    inline AcdGapId()
      :m_val(0){
    }
    
    inline AcdGapId(int type, int gap, int face, int row, int col)
      :m_val(0){
      setVal(type,gap,face,row,col);
    }

    inline AcdGapId(unsigned int val)
      :m_val(val){;}
    
    inline AcdGapId(const AcdGapId& other)
      :m_val(other.m_val){;}
    
    inline ~AcdGapId(){;}
    
    inline AcdGapId& operator=(const AcdGapId& other) {
      m_val = other.asShort();
      return *this;
    }
    inline bool operator==(const AcdGapId& other) {
      return m_val == other.asShort();
    }  
    inline unsigned short asShort() const {
      return m_val;
    }
    inline unsigned char gapType() const {
      return ( NibbleMask & ( m_val >> TypeShift ) );
    }
    inline unsigned char face() const {
      return ( ThreeBitMask & ( m_val >> FaceShift ) );
    }
    inline unsigned char row() const {
      return ( ThreeBitMask & ( m_val >> RowShift ) );
    }    
    inline unsigned char col() const {
      return ( ThreeBitMask & ( m_val >> ColShift ) );
    }
    inline unsigned char gap() const {
      return ( ThreeBitMask & ( m_val >> GapShift ) );
    }    
    inline unsigned short closestTile() const {
      return 100*face() + 10*row() + col();
    }
    inline unsigned int asDecimal() const {
      return 10000*gapType() + 1000*gap() + 100*face() + 10*row() + col();
    }
    void setVal(unsigned char type, unsigned char gap, unsigned char face, unsigned char row, unsigned char col) {
      m_val = 0;
      m_val |= ( NibbleMask & type ) << TypeShift;
      m_val |= ( ThreeBitMask & gap ) << GapShift;
      m_val |= ( ThreeBitMask & face ) << FaceShift;
      m_val |= ( ThreeBitMask & row ) << RowShift;
      m_val |= ( ThreeBitMask & col ) << ColShift;  
    }

  private:
    
    unsigned short int m_val;
    
  };

};

#endif
