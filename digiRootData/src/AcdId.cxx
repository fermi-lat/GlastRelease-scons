///////////////////////////////////////////////////////////////////////////
//									  
// The AcdId class contains a tile ID number.	
// AcdId is used in the AcdTile class.
//									  
///////////////////////////////////////////////////////////////////////////

// Apr 2000 Daniel Flath - Minor changes to function names, etc.
// Jan 1999 Daniel Flath - ROOT HTML Documentation added
// Dec 1999 Daniel Flath - Rewrite for GLAST
// Version 1.1 25 Oct 1999 R.Dubois Clone from LCD towerID

#include "digiRootData/AcdId.h"

ClassImp(AcdId)

UShort_t AcdId::badId = 3799;

///________________________________________________________________________
AcdId::AcdId() : m_id(0), m_used(1) {
}

AcdId::AcdId(const AcdId& id) : m_id(id.m_id), m_used(id.m_used) 
{ 
};

AcdId::AcdId(short l, short f, short r, short c) {
    m_used = 1;
    m_id = 0;
    setLayer(l);
    setFace(f);
    setRow(r);
    setColumn(c);
};

AcdId::AcdId(UInt_t i, short base, short used) {
    setId(i, used, base);
}

///________________________________________________________________________
UInt_t AcdId::getId(short base) const { 
    // Returns the tile ID number
    if (base == 2) return m_id; 
    return (getLayer() * 1000 + getFace() * 100 + getRow() * 10 + getColumn());
}

///________________________________________________________________________
void AcdId::setId(UInt_t newVal, short used, short base) { 
// Sets the ACD ID number to newVal
// If this is an unused PMT, store that information
    if (used <= 0) {
        m_used = 0;
        m_id = 0;
        setLayer(3);
        setFace(7);
        setRow(9);
        setColumn(9);
    } else {
        m_used = 1;
        if (base == 2) {
            m_id = newVal;
        } else if (base == 10) {
            short lay, face, row, col;
            base10ToAcdId(newVal, lay, face, row, col);
            setLayer(lay);
            setFace(face);
            setRow(row);
            setColumn(col);
        }
    }
}


bool AcdId::isTop () const 
{ 
    return (getFace() == 0); 
}

bool AcdId::isSide () const 
{ 
    return (getFace() != 0); 
}

short AcdId::getLayer () const 
{ 
    return ((m_id & _layermask) >> layerShift); 
}

short AcdId::getFace () const
{ 
    return (word(2, (m_id & _facemask))); 
}

short AcdId::getRow () const 
{ 
    return(word (1, m_id)); 
}

short AcdId::getColumn () const 
{ 
    return(word (0, m_id)); 
}

void AcdId::setLayer( unsigned int val)
{ 
    m_id &= ~_layermask;
    m_id |= (val << layerShift);
    //short two = 2;
    
   // set_word( two, m_id, 
   //     ((val == 0) ? 0 : 8) | ((_facemask & m_id) >> 8) );
}

void AcdId::setFace( unsigned int f )
{
    short two = 2;
    set_word( two, m_id, 
        f | ((_layermask & m_id ) >> 8) );
}

void AcdId::setRow( unsigned int r ) 
{ 
    short one = 1;
    set_word( one, m_id, r );
}

void AcdId::setColumn( unsigned int c ) 
{ 
    short zero = 0;
    set_word( zero, m_id, c ); 
}

bool AcdId::wasConnected() const
{
    return ( (m_used == 1) ? true : false);
}

void AcdId::setConnected(Char_t c) {
    m_used = c;
}

void AcdId::base10ToAcdId(UInt_t val, short &lay, short &face, short &row, short &col) {
    lay = val / 1000;
    val -= lay*1000;
    face = val / 100;
    val -= face*100;
    row = val / 10;
    val -= row*10;
    col = val;
}