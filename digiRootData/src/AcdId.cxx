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


///________________________________________________________________________
UShort_t AcdId::getId() const { 
    // Returns the tile ID number
    
    return m_id;
}
///________________________________________________________________________
void AcdId::setId(Short_t newVal) { 
// Sets the ACD ID number to newVal
// If this is an unused PMT, store that information
    if (newVal < 0) {
        m_used = 0;
        m_id = 0;
    } else {
        m_used = 1;
        m_id = UShort_t(newVal);
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
    return ((m_id & _layermask) ? 1 : 0); 
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
    short two = 2;
    
    set_word( two, m_id, 
        ((val == 0) ? 0 : 8) | ((_facemask & m_id) >> 8) );
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