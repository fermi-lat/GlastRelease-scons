
#ifndef AcdId_H
#define AcdId_H

#include "TObject.h"

/*! \class AcdId
\brief encapsulate the id for an ACD tile.

 ACD tile numbering
 Layer 0 or 1, where 0 corresponds to the inner layer
 Face 0 - 4 
          Face 0 == top  (hat)
          Face 1 == -X side 
          Face 2 == -Y side
          Face 3 == +X side
          Face 4 == +Y side
 Row 0 - 9
          On Face 0 (top) rows are numbering along +x (or +y) axis
          On the sides, the rows are numbered from the front of the instrument (closest to the top)

 Column 0 - 9
          On the sides, columns are numbered along +x (or +y) axis

 One could imagine that an ACD id could be represented by 12 bits.
 
  __       __ __ __     __ __ __ __      __ __ __ __
 LAYER       FACE           ROW            COLUMN
 where the layer is the most significant bit


Jun 2001 Heather Kelly - renamed from TileID to AcdId
Apr 2000 Daniel Flath - Minor changes to function names, etc.
Jan 1999 Daniel Flath - ROOT HTML Documentation added
Dec 1999 Daniel Flath - Rewrite for GLAST
Version 1.1 25 Oct 1999 R.Dubois Clone from LCD towerID
*/

class AcdId : public TObject {

private:
    /// Proper Ritz index number for this tile
    UInt_t m_id;  
    /// Denotes whether the PMT was truly connected to a tile
    Char_t m_used;
    /// Readout order - starting from 0
    UShort_t m_readoutIndex;
    /// set layer
    void setLayer( unsigned int val );
    /// set the face number
    void setFace (unsigned int f);
    /// set the row
    void setRow( unsigned int r );
    /// set the column
    void setColumn( unsigned int c );
    
    enum {
        _layermask = 0x0800,
        _facemask  = 0x0700,
        _rowmask   = 0x00F0,
        _colmask   = 0x000F
    };
    
public:
    AcdId();
    AcdId(const AcdId& id);
    AcdId(short l, short f, short r, short c);
    virtual ~AcdId() { };
    
    UInt_t getId() const;
    void setId(Short_t newVal);
    /// is this a top tile?
    bool isTop () const;  
    /// is this a side tile?
    bool isSide () const;   
    /// which layer is this tile in (1st(0 - closer to center of GLAST) or 2nd (1))
    short getLayer () const;     
    /// which face number?
    short getFace() const;
    /// which row number?
    short getRow() const;     
    /// which column?
    short getColumn () const;   
    /// was this PMT actually connected to a tile?
    bool wasConnected() const;
    /// set the bit to denote whether this PMT was connected 0 == not connected,1 is connected
    void setConnected(Char_t c);
    /// retrieve the readout index
    UShort_t getReadoutIndex() const { return m_readoutIndex; };
    /// set the readout index for this tile
    void setReadoutIndex(UShort_t i) { m_readoutIndex = i; };
    
    static UShort_t badId;
    
private:
    inline short int word ( short i, const UInt_t& v ) const
    {
        short   shft = i*4;
        unsigned int  mask = 0x000F << shft;
        return (v & mask) >> shft;
    }
    
    void set_word ( short i, UInt_t& v, const UInt_t& s )
    {
        short   shft = i*4;
        unsigned int  mask = 0x000F << shft;
        v = ((v & ~mask) | ((s << shft) & mask));
    }
    
    ClassDef(AcdId,3)	  // ACD Tile identification information
};

#endif
