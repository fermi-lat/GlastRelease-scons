// $Id$
//
// Interface for the ScintillatorId class.

#ifndef _H_Scintillator_Id
#define _H_Scintillator_Id

#include "idents/ModuleId.h"

/** ScintillatorId - class to encapsulate the id for a veto-tile.

 From Heather Arrighi:
 ---------------------
 Id numbers for the Tiles look like this:
 2LayerBit Top/SideBit NSEW(2bits) HiOrderId(4 bits) LoOrderId(4 bits)
 2LayerBit == 0 is first layer, 1 is 2nd layer
 Top/SideBit == 0 for top tile, 1 for side tile
 EWNS == 1st bit denotes East(0) or West(1) Side       NOTE only used for Side tiles
         2nd bit denotes North(0) or South(1) Side        will be 0 for Top tiles
 Top Tiles:   HiOrderId == position along E/W direction  W == 0000
              LoOrderId == position along N/S direction  S == 0000
 Side Tiles:  HiOrderId == Row Number where 0 == top row of tiles
              LoOrderId == tile Number (0 - ?) From S to N OR W to E

*/
namespace idents {
class   ScintillatorId {
public:
    ScintillatorId ();
    ScintillatorId ( const ScintillatorId& );
    ScintillatorId ( unsigned int id ) : m_id (id) {}
    virtual ~ScintillatorId ();

    // access
    operator const unsigned int& () const { return m_id; }

    void write (std::ostream& out) const
    {
        out << m_id;
    }

    void read (std::istream& in)
    {
        in >> m_id;
    }

    // all tile access methods
    inline bool    top () const;       // is this a top tile?
    inline bool    side () const;      // is this a side tile?
    inline short   layer () const;     // which layer is this tile in (0-1)

    // top tile access methods
    inline short   row () const;       // N-S row (N = 0)
    inline short   column () const;    // E-W column (W = 0)
    
    // side tile access methods
    inline short   vertical () const;      // vertical position (0 top)
    inline short   orientation () const;   // NWSE orientation (N = 0)
    inline short   position () const;      // N-S or E-W position on side (N,W = 0)

    // is this tile assoicated to a particular tower module ?
    // returns a bit pattern with the bit set representing the position of the ACD tile as follows:
    // [4 bits = top] [4 bits = W] [4 bits = S] [4 bits = E] [4 bits = N] 
    // where the bit that is set within a particular group of 4 corresponds to the position
    // of the tile starting in the upper left and moving clockwise while looking either straight down,
    // due north or due west through the instrument.
    int association ( const idents::ModuleId& ) const;

protected:

    enum {xNum=4, yNum=4}; // assumed dimension of towers
    // all tiles
    inline void top (bool t);           // set to top tile
    inline void side (bool t);          // set to side tile
    inline void layer (short l);        // set layer number

    // top tile settings
    inline void row (short r);      // set N-S row (N = 0)
    inline void column (short c);   // set E-W column (W = 0)
    
    // side tile access methods
    inline void vertical (short v);     // set vertical pos on side (0 = bottom)
    inline void orientation (short o);  // set NWSE orientation (N = 0)
    inline void position (short p);     // set N-S/E-W position in row (N,W = 0)

    // friends
    friend class ACDTopMed;
    friend class ACDSideMed;
    friend class Scintillator;

private:
    // private access methods
    enum {
        _layermask = 0x0800,
        _sidemask = 0x0400,
        _neswmask = 0x0300, 
    };

    unsigned int    m_id;   // id value (4 byte word)
};

// inline declarations
namespace bitmanip {

// template functions for getting/setting four bit words
// i is the word #
// v is the value
// s it the new value

template <class _T>
inline short int   word ( short i, const _T& v )
{
    short   shft = i*4;
    _T  mask = 0x000F << shft;
    return (v & mask) >> shft;
}

template <class _T>
void    set_word ( short i, _T& v, const _T& s )
{
    short   shft = i*4;
    _T  mask = 0x000F << shft;
    v = ((v & ~mask) | ((s << shft) & mask));
}

} // namespace bitmanip


inline ScintillatorId::ScintillatorId () : m_id (0) {}
inline ScintillatorId::ScintillatorId ( const ScintillatorId& r ) : m_id (r.m_id) {}
inline ScintillatorId::~ScintillatorId () {}

inline bool ScintillatorId::top () const 
{ return !side(); }

inline bool ScintillatorId::side () const 
{ return ((m_id & _sidemask) != 0x0000); }

inline short ScintillatorId::layer () const 
{ return (m_id & _layermask) ? 1 : 0; }

inline short ScintillatorId::row () const 
{ return bitmanip::word (1, m_id); }

inline short ScintillatorId::column () const 
{ return bitmanip::word (0, m_id); }

inline short ScintillatorId::vertical () const 
{ return bitmanip::word (1, m_id); }

inline short ScintillatorId::orientation () const 
{ return ( bitmanip::word (2, m_id & _neswmask) ); }

inline short ScintillatorId::position () const 
{ return bitmanip::word (0, m_id); }

inline void ScintillatorId::top (bool t)
{ m_id = ((~_sidemask & m_id) | ((!t) ? _sidemask : 0)); }

inline void ScintillatorId::side (bool t)
{ top (!t); }

inline void ScintillatorId::layer (short l)
{ bitmanip::set_word<unsigned int> (2, m_id, l | bitmanip::word (2, ~_layermask & m_id));  }

inline void ScintillatorId::row (short r) 
{ bitmanip::set_word<unsigned int> (1, m_id, r); }

inline void ScintillatorId::column (short c) 
{ bitmanip::set_word<unsigned int> (0, m_id, c); }

inline void ScintillatorId::vertical (short v) 
{ bitmanip::set_word<unsigned int> (1, m_id, v); }

inline void ScintillatorId::orientation (short o) 
{ bitmanip::set_word<unsigned int> (2, m_id, o | bitmanip::word (2, ~_neswmask & m_id)); }

inline void ScintillatorId::position (short p)
{ bitmanip::set_word<unsigned int> (0, m_id, p); }

inline int  ScintillatorId::association ( const idents::ModuleId& m ) const
{
    // decide whether or not the veto tile is
    // in proximity to the tower
    int rc = 0;     // return codes
    int lpos = 0;   // 4-bit value for the specific tile position
    unsigned int shft = 0;  // how many bits to shift the 4-bit group
    int vy, vx, ty = m.iy() - 1, tx = m.ix() - 1;
    if (side()) {
        vy = orientation();
        tx = (orientation() & 1) ? m.iy() - 1 : m.ix() - 1;
        vx = position();

        // number of bits to slide to the left to assign the proper 4-bit group
        shft = vy * 4;
                
        // check for min/max in row v tile orientation
        // this makes sure the tower is oriented on the appropriate edge of the 
        // instrument.
        if ((vy == 2) || (vy == 1)) {
            if (vy & 1) {
                if (m.ix() > 1)   return 0;
            } else {
                if (m.iy() > 1)   return 0;
            }
        } else {
            if (vy & 1) {
                if (m.ix() != xNum) return 0;
            } else {
                if (m.iy() != yNum) return 0;
            }
        }

        // compare row/column indexes
        if (vx == tx) {
            // if the ACD is tower-aligned - set all 4 bits high for this word.
            //if (Tower::align)   lpos = 0x000F;  
            // otherwise this is postion 0 or 2

        } else if (/*!Tower::align &&*/ (vx == (tx + 1))) {
            // this is position 1 or 3
            lpos |= 1;

        } else
            return 0;

        // add 2 to the upper row of ACD tiles on side
        lpos |= (vertical())*2;
        // shift the bits to the proper 4 bit position

    } else {
        // top tile - highest position...
        shft = 16;

        // figure out the x-position of the tile
        vx = column();
        if (vx == tx) { // is the acd column == the tower column?
            // if the ACD is tower-aligned - set all 4 bits high for this word.
            //if (Tower::align)   lpos = 0x000F; 
            // otherwise this is postion 0 or 2

        } else if (/*!Tower::align &&*/ (vx == (tx + 1))) {
            // this is position 1 or 3
            lpos |= 1;
        } else  
            return 0;   // not associated
        
        // figure out the y position of the tile
        vy = row();
        if (vy == ty) {
            // if the ACD is tower-aligned - set all 4 bits high for this word.
            //if (Tower::align)   lpos = 0x000F;
            // this is position  0 or 1

        } else if (/*!Tower::align && */(vy == (ty + 1))) {
            // this is postition 2 or 3
            lpos |= 2;

        } else
            return 0;   // not associated

    }
    // slide the bits to the left by 16 to the top position
    rc = 1 << lpos;
    rc = rc << shft;
    return rc;
}

inline std::ostream& operator<<(std::ostream& out,const ScintillatorId& id){id.write(out); return out;}
inline std::istream& operator>>(std::istream&  in, ScintillatorId& id){id.read(in); return in;}
} // namespace idents
#endif // _H_Scintillator_Id
