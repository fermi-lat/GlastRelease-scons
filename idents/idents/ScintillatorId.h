// $Id$
//
// Interface for the ScintillatorId class.

#ifndef _H_Scintillator_Id
#define _H_Scintillator_Id

#include "facilities/bitmanip.h"
#include "idents/ModuleId.h"
#include "instrument/Tower.h"
#include "instrument/Glast.h"

#include <vector>

/*!
\class ScintillatorId 
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
*/

namespace idents {
class   ScintillatorId {
public:
    ScintillatorId ();
    ScintillatorId ( const ScintillatorId& );
    ScintillatorId ( unsigned int id ) : m_id (id) { m_towers.clear(); }
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
    inline unsigned int id ();
    inline bool    top () const;       // is this a top tile?
    inline bool    side () const;      // is this a side tile?
    inline short   layer () const;     // which layer is this tile in (0-1)

    inline short   face() const;
    inline short   row() const;     
    inline short   column () const;   
    std::vector<idents::ModuleId> getAssociatedTowers() { return m_towers; }

    /// add a new ModId to the list of associated towers
    inline void addTower(const idents::ModuleId &m); 

    /// is this tile assoicated to a particular tower module ?
    /// returns a bit pattern with the bit set representing the position of the ACD tile as follows:
    /// [4 bits = top] [4 bits = W] [4 bits = S] [4 bits = E] [4 bits = N] 
    /// where the bit that is set within a particular group of 4 corresponds to the position
    /// of the tile starting in the upper left and moving clockwise while looking either straight down,
    /// due north or due west through the instrument.
    int association ( const idents::ModuleId &m ) const;

    /// determines whether two tile are adjacent
    bool neighbor( const ScintillatorId& ) const;

protected:

    inline void layer( unsigned int val ); // set layer
    inline void face (unsigned int f);   // set face 
    inline void row( unsigned int r );   // set row 
    inline void column( unsigned int c ); // set column
    
    // friends
    friend class ACDTopMed;
    friend class ACDSideMed;
    friend class Scintillator;

private:

    enum {
        _layermask = 0x0800,
        _facemask  = 0x0700,
        _rowmask   = 0x00F0,
        _colmask   = 0x000F
    };

    unsigned int    m_id;   // id value (4 byte word)

    // list of modules Ids associated with this ScintillatorID
    std::vector<idents::ModuleId> m_towers;
};

// inline declarations

inline ScintillatorId::ScintillatorId () : m_id (0) {m_towers.clear();}
inline ScintillatorId::ScintillatorId ( const ScintillatorId& r ) : m_id (r.m_id), m_towers(r.m_towers) {}
inline ScintillatorId::~ScintillatorId () {m_towers.clear();}

inline unsigned int ScintillatorId::id() 
{ return (layer() * 1000 + face() * 100 + row() * 10 + column()); }

inline bool ScintillatorId::top () const 
{ return (face() == 0); }

inline bool ScintillatorId::side () const 
{ return (face() != 0); }

inline short ScintillatorId::layer () const 
{ return (m_id & _layermask) ? 1 : 0; }

inline short ScintillatorId::face () const
{ return bitmanip::word(2, (m_id & _facemask)); }

inline short ScintillatorId::row () const 
{ return bitmanip::word (1, m_id); }

inline short ScintillatorId::column () const 
{ return bitmanip::word (0, m_id); }

inline void ScintillatorId::layer( unsigned int val)
{ 
  short two = 2;

  bitmanip::set_word<unsigned int>( 
      two, m_id, 
      ((val == 0) ? 0 : 8) | ((_facemask & m_id) >> 8) );
}

inline void ScintillatorId::face( unsigned int f )
{
    short two = 2;
    bitmanip::set_word<unsigned int> (
        two, m_id, 
        f | ((_layermask & m_id ) >> 8) );
}

inline void ScintillatorId::row( unsigned int r ) 
{ 
  short one = 1;
  bitmanip::set_word<unsigned int>( one, m_id, r );
}

inline void ScintillatorId::column( unsigned int c ) 
{ 
  short zero = 0;
  bitmanip::set_word<unsigned int>( zero, m_id, c ); 
}

inline void ScintillatorId::addTower(const idents::ModuleId &m) {
    m_towers.push_back(m);
}

inline int  ScintillatorId::association ( const idents::ModuleId& m ) const
{
    typedef std::vector < idents::ModuleId >::const_iterator const_iterator;

    for( const_iterator it = m_towers.begin(); it != m_towers.end(); ++it) {
        const idents::ModuleId& modId = *it;
	if (modId == m) return 1;
    }
    
    return 0;
}

inline bool ScintillatorId::neighbor( const ScintillatorId& id) const
{
    // tiles are on the same face
    if (id.face() == face()) {
        if ( (fabs(id.row() - row()) <= 1) && (fabs(id.column() - column()) <= 1) )
            return true; 
    } else {
        // otherwise tiles are on different faces
       // if (column() == 0) || (column() == 4) {

        //}

    }

    return false;
}

inline std::ostream& operator<<(std::ostream& out,const ScintillatorId& id){id.write(out); return out;}
inline std::istream& operator>>(std::istream&  in, ScintillatorId& id){id.read(in); return in;}
} // namespace idents

#endif // _H_Scintillator_Id
