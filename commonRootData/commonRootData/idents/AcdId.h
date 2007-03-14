
#ifndef commonRootData_AcdId_H
#define commonRootData_AcdId_H

#include "TObject.h"
#include <string>

/** @class AcdId
 @brief encapsulate the id for an ACD detector.
 
 @verbatim
 ACD tile numbering
 N/A 0 - 1, where 1 corresponds to a non-connected PMT
 Face 0 - 4 
    Face 0 == top  (hat)
    Face 1 == -X side 
    Face 2 == -Y side
    Face 3 == +X side
    Face 4 == +Y side
 Row 0 - 9
    On Face 0 (top) rows are numbering along +x (or +y) axis
    On the sides, the rows are numbered from the front of the instrument 
    (closest to the top)
 Column 0 - 9
    On the sides, columns are numbered along +x (or +y) axis

 One could imagine that an ACD id could be represented by 13 bits.
 
 __  __       __ __ __     __ __ __ __      __ __ __ __
 N/A         FACE           ROW            COLUMN
              OR                              OR
            RIBBON ORIENTATION              RIBBON NUMBER
 where the layer bits are the most significant bits
@endverbatim

Some routines in this class have a default parameter denoted as "base"
This refers to those instances where an input parameter could be in 
base 10 or base 2.  In all cases, the default is base 10.  If a user
desires to use base 2, the parameter called "base" should be supplied
and set equal to 2.

*/

class AcdId : public TObject {

public:

  //
  static void convertToTilePmt(UInt_t  cable, UInt_t  channel, UInt_t& tile, UInt_t& pmt);
  static void convertToGarcGafe(UInt_t tile, UInt_t pmt, UInt_t& garc, UInt_t& gafe);
  static UInt_t           tileFromGemIndex(UInt_t gemIndex);
  static UInt_t           gemIndexFromTile(UInt_t tile);
  static UInt_t           tileFromIndex(UInt_t gemIndex);
  static UInt_t           indexFromTile(UInt_t tile);
  
private:
    // extend the notion of a layer to allow for more than 2 layers
    // for XGTs for example
    enum {
        _layermask = 0x1800,
        _namask = 0x1800,
        _facemask  = 0x0700,
        _rowmask   = 0x00F0,
        _colmask   = 0x000F,
        layerShift = 11,
        naShift = 11,
        _ribbonmask = 0x0007,
        _ribbonorientmask = 0x0700,
        maxAcdTileFace = 4,
        ribbonX = 5,  // ribbons that extend along x-axis
        ribbonY = 6   // ribbons that extend along y-axis
    };
    
public:

    AcdId();
    AcdId(const AcdId& id);
    AcdId(UShort_t l, UShort_t f, UShort_t r, UShort_t c);
    AcdId(UInt_t i, Short_t base=10, Short_t used=1);
    AcdId(UShort_t ribOriet, UShort_t ribNum);
    
    void initialize(UShort_t l, UShort_t f, UShort_t r, UShort_t c);
    void initialize(UShort_t ribOrient, UShort_t ribNum);

    virtual ~AcdId() { };

    void Clear(Option_t *option ="");
    Bool_t CompareInRange( const AcdId &, const std::string & ="unknown AcdId" ) const ; // for tests
    void Print(Option_t *option="") const;

    /// retrieve the id, either in base 2
    /// or in base 10 denoting the matrix number
    UInt_t getId(Short_t base = 10) const;
    /// set the id, specifying whether this is in base 10 or base 2
    void setId(UInt_t newVal, Short_t used = 1, Short_t base=10);

    /// is this a tile?
    bool isTile() const;
    /// is this a ribbon?
    bool isRibbon() const;
    /// is this a top tile?
    bool isTop () const;  
    /// is this a side tile?
    bool isSide () const;   
    /// which layer is this tile in (1st(0 - closer to center of GLAST) or 2nd (1))
    Short_t getLayer () const;     
    /// Is this not connected?
    Short_t getNa () const;     
    /// which face number?
    Short_t getFace() const;
    /// which row number?
    Short_t getRow() const;     
    /// which column?
    Short_t getColumn () const;   
    /// Number of the ribbon starting from 0 and increasing along increasing X or Y
    Short_t getRibbonNumber() const;
    /// Does this ribbon extend along the X or Y axis?
    Short_t getRibbonOrientation() const;

    /// was this PMT actually connected to a tile?
    bool wasConnected() const;
    /// set the bit to denote whether this PMT was connected 0 == not connected,1 is connected
    void setConnected(Char_t c);
    /// retrieve the readout index
    UShort_t getReadoutIndex() const { return m_readoutIndex; };
    /// set the readout index for this tile
    void setReadoutIndex(UShort_t i) { m_readoutIndex = i; };
    
    /// convert a value in base 10, to the appropriate layer, face, row, column
    void base10ToAcdId(UInt_t val, Short_t &lay, Short_t &face, 
        Short_t &row, Short_t &col);

	bool operator==(const AcdId& i) const;
    
    /// Special id value to denote an improper id number
    static UShort_t badId;
    
private:

    inline short int word ( Short_t i, const UInt_t& v ) const
    {
        Short_t   shft = i*4;
        UInt_t  mask = 0x000F << shft;
        return (v & mask) >> shft;
    }
    
    void set_word ( Short_t i, UInt_t& v, const UInt_t& s )
    {
        Short_t   shft = i*4;
        UInt_t  mask = 0x000F << shft;
        v = ((v & ~mask) | ((s << shft) & mask));
    }

    /// set layer
    void setLayer(UInt_t val);
    /// set N/A
    void setNa(UInt_t val);
    /// set the face number
    void setFace (UInt_t f);
    /// set the row
    void setRow(UInt_t r);
    /// set the column
    void setColumn(UInt_t c);
    
    void setRibbonNumber(UInt_t r);
    void setRibbonOrientation(UInt_t r);

    /// Proper Ritz index number for this tile
    UInt_t m_id;  
    /// Denotes whether the PMT was truly connected to a tile
    Char_t m_used;
    /// Readout order - starting from 0
    UShort_t m_readoutIndex;
    
    ClassDef(AcdId,4)	  // ACD Tile identification information
};

#endif
