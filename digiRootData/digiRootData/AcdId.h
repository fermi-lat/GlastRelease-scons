
#ifndef AcdId_H
#define AcdId_H

#include "TObject.h"

/*! \class AcdId
\brief encapsulate the id for an ACD tile.

\verbatim
 ACD tile numbering
 Layer 0 - 3, where 0 corresponds to the inner layer
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

 One could imagine that an ACD id could be represented by 13 bits.
 
 __  __       __ __ __     __ __ __ __      __ __ __ __
 LAYER       FACE           ROW            COLUMN
 where the layer bits are the most significant bits
\endverbatim

Some routines in this class have a default parameter denoted as "base"
This refers to those instances where an input parameter could be in 
base 10 or base 2.  In all cases, the default is base 10.  If a user
desires to use base 2, the parameter called "base" should be supplied
and set equal to 2.

\li Jun 2001 Heather Kelly - renamed from TileID to AcdId
\li Apr 2000 Daniel Flath - Minor changes to function names, etc.
\li Jan 1999 Daniel Flath - ROOT HTML Documentation added
\li Dec 1999 Daniel Flath - Rewrite for GLAST
\li Version 1.1 25 Oct 1999 R.Dubois Clone from LCD towerID
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
    
    // extend the notion of a layer to allow for more than 2 layers - for XGTs for example
    enum {
        _layermask = 0x1800,
        _facemask  = 0x0700,
        _rowmask   = 0x00F0,
        _colmask   = 0x000F,
        layerShift = 11
    };
    
public:
    AcdId();
    AcdId(const AcdId& id);
    AcdId(short l, short f, short r, short c);
    AcdId(UInt_t i, short base=10, short used=1);
    virtual ~AcdId() { };
    
    /// retrieve the id, either in base 2, or in base 10, denoting the matrix number
    UInt_t getId(short base = 10) const;
    /// set the id, specifying whether this is in base 10 or base 2
    void setId(UInt_t newVal, short used = 1, short base=10);
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
    /// convert a value in base 10, to the appropriate layer, face, row, column
    void base10ToAcdId(UInt_t val, short &lay, short &face, short &row, short &col);
    
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
