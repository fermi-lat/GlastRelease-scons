//									  
// The AcdId class contains a tile ID number.	
// AcdId is used in the AcdTile class.
//									  

// Apr 2000 Daniel Flath - Minor changes to function names, etc.
// Jan 1999 Daniel Flath - ROOT HTML Documentation added
// Dec 1999 Daniel Flath - Rewrite for GLAST
// Version 1.1 25 Oct 1999 R.Dubois Clone from LCD towerID

#include "digiRootData/AcdId.h"
#include <iostream>

ClassImp(AcdId)

UShort_t AcdId::badId = 3799;

AcdId::AcdId() : m_id(0), m_used(1) {
}

AcdId::AcdId(const AcdId& id) : m_id(id.m_id), m_used(id.m_used) 
{ 
}

AcdId::AcdId(UShort_t l, UShort_t f, UShort_t r, UShort_t c) {
    m_used = 1;
    m_id = 0;
    setLayer(l);
    setFace(f);
    setRow(r);
    setColumn(c);
}

AcdId::AcdId(UInt_t i, Short_t base, Short_t used) {
    setId(i, used, base);
}

void AcdId::Clear(Option_t *option) {
    m_used = 1;
    m_id = 0;
}

void AcdId::Print(Option_t *option) const {
    using namespace std;
    TObject::Print(option);
    cout.precision(2);
    cout << "Id: " << getId()
        << "  (Layer, Face, Row, Column): (" << getLayer() << ", "
        << getFace() << ", " << getRow() << ", " << getColumn()
        << ")" << endl;
}

UInt_t AcdId::getId(Short_t base) const { 
    // Returns the tile ID number
    if (base == 2) return m_id; 
    return (getLayer() * 1000 + getFace() * 100 + getRow() * 10 + getColumn());
}

void AcdId::setId(UInt_t newVal, Short_t used, Short_t base) { 
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

Short_t AcdId::getLayer () const 
{ 
    return ((m_id & _layermask) >> layerShift); 
}

Short_t AcdId::getFace () const
{ 
    return (word(2, (m_id & _facemask))); 
}

Short_t AcdId::getRow () const 
{ 
    return(word (1, m_id)); 
}

Short_t AcdId::getColumn () const 
{ 
    return(word (0, m_id)); 
}

void AcdId::setLayer( UInt_t val )
{ 
    m_id &= ~_layermask;
    m_id |= (val << layerShift);
    //short two = 2;
    
   // set_word( two, m_id, 
   //     ((val == 0) ? 0 : 8) | ((_facemask & m_id) >> 8) );
}

void AcdId::setFace( UInt_t f )
{
    short two = 2;
    set_word( two, m_id, 
        f | ((_layermask & m_id ) >> 8) );
}

void AcdId::setRow( UInt_t r ) 
{ 
    short one = 1;
    set_word( one, m_id, r );
}

void AcdId::setColumn( UInt_t c ) 
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

void AcdId::base10ToAcdId(UInt_t val, Short_t &lay, Short_t &face, 
                          Short_t &row, Short_t &col) {
    lay = val / 1000;
    val -= lay*1000;
    face = val / 100;
    val -= face*100;
    row = val / 10;
    val -= row*10;
    col = val;
}

bool AcdId::operator== (const AcdId& i) const {
	if (m_id == i.m_id) return true;
	return false;
}