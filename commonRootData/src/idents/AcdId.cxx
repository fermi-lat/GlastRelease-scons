//									  
// The AcdId class contains a tile ID number.	
// AcdId is used in the AcdTile class.
//									  

// Apr 2000 Daniel Flath - Minor changes to function names, etc.
// Jan 1999 Daniel Flath - ROOT HTML Documentation added
// Dec 1999 Daniel Flath - Rewrite for GLAST
// Version 1.1 25 Oct 1999 R.Dubois Clone from LCD towerID

#include <commonRootData/idents/AcdId.h>
#include <commonRootData/RootDataUtil.h>
#include <iostream>

ClassImp(AcdId)

UShort_t AcdId::badId = 3799;

AcdId::AcdId()
 : m_id(0), m_used(1) {
}

AcdId::AcdId(const AcdId& id)
 : TObject(id), m_id(id.m_id), m_used(id.m_used) 
{ 
}

AcdId::AcdId(UShort_t l, UShort_t f, UShort_t r, UShort_t c)
 : m_id(0), m_used(1) {
    m_used = 1;
    m_id = 0;
    setLayer(l);
    setFace(f);
    setRow(r);
    setColumn(c);
}

AcdId::AcdId(UShort_t ribbonOrient, UShort_t r)
 : m_id(0), m_used(1) {
    m_used = 1;
    m_id = 0;
    setRibbonNumber(r);
    setRibbonOrientation(ribbonOrient);
}


AcdId::AcdId(UInt_t i, Short_t base, Short_t used)
 : m_id(0), m_used(1) {
    setId(i, used, base);
}

void AcdId::initialize(UShort_t l, UShort_t f, UShort_t r, UShort_t c) {
    m_used = 1;
    m_id = 0;
    setLayer(l);
    setFace(f);
    setRow(r);
    setColumn(c);
}

void AcdId::initialize(UShort_t ribOrient, UShort_t ribNum) {
    m_used = 1;
    m_id = 0;
    setRibbonNumber(ribNum);
    setRibbonOrientation(ribOrient);
}

void AcdId::Clear(Option_t *) {
    m_used = 1;
    m_id = 0;
}

Bool_t AcdId::CompareInRange( const AcdId & id, const std::string & name ) const {

    bool result = true ;

    result = rootdatautil::CompareInRange(m_id,id.m_id,"Id") && result ;
    result = rootdatautil::CompareInRange(m_used,id.m_used,"Used") && result ;
    // apparently readoutIndex is obsolete
    //result = rootdatautil::CompareInRange(m_readoutIndex,id.m_readoutIndex,"ReadoutIndex") && result ;

    if (!result) {
        std::cout<<"Comparison ERROR for "<<name<<std::endl ;
    }
    return result ;

}

void AcdId::Print(Option_t *option) const {
    using namespace std;
    TObject::Print(option);
    cout.precision(2);
    cout << "Id: " << getId() << endl;
    if (isTile())
        cout << "  ( Face, Row, Column): (" 
        << getFace() << ", " << getRow() << ", " << getColumn()
        << ")" << endl;
    else if (isRibbon())
        cout << "   (RibbonNum, RibbonOrient): (" << getRibbonNumber() << ", "
        << getRibbonOrientation() << ") " << endl;
    else
        cout << " N/A: " << getNa() << endl;
}

UInt_t AcdId::getId(Short_t base) const { 
    // Returns the tile ID number
    if (base == 2) return m_id; 
    if (isTile())
        return (getNa() * 1000 + getFace() * 100 + getRow() * 10 + getColumn());
    else 
        return (getNa() * 1000 + getRibbonOrientation() * 100 + getRibbonNumber());
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

bool AcdId::isTile () const
{ return ((!getNa()) && getFace() <= maxAcdTileFace); }

bool AcdId::isRibbon () const
{ return ((!getNa()) && getFace() > maxAcdTileFace); }

bool AcdId::isTop () const 
{ 
    return ((!getNa()) && getFace() == 0); 
}

bool AcdId::isSide () const 
{ 
    return ((!getNa()) && !isRibbon() && (getFace() != 0)); 
}

Short_t AcdId::getLayer () const 
{ 
    return (getNa()); 
}

Short_t AcdId::getNa () const 
{ 
    return ((m_id & _namask) >> naShift); 
}

Short_t AcdId::getFace () const
{ 
    return (word(2, (m_id & _facemask))); 
}

Short_t AcdId::getRow () const 
{ 
    if (isTile())
        return(word (1, m_id)); 
    return -1;
}

Short_t AcdId::getColumn () const 
{ 
    if (isTile())
        return(word (0, m_id)); 
    return -1;
}

Short_t AcdId::getRibbonNumber () const
{ 
    if (isRibbon()) return word(0, (m_id & _ribbonmask)); 
    return -1;
}

Short_t AcdId::getRibbonOrientation () const
{
    if (isRibbon()) return word(2, (m_id & _ribbonorientmask));
    return -1;
}


void AcdId::setLayer( UInt_t val )
{ 
    setNa(val);
}

void AcdId::setNa( UInt_t val )
{ 
    m_id &= ~_namask;
    m_id |= (val << naShift);
    //short two = 2;
    
   // set_word( two, m_id, 
   //     ((val == 0) ? 0 : 8) | ((_facemask & m_id) >> 8) );
}

void AcdId::setFace( UInt_t f )
{
    short two = 2;
    set_word( two, m_id, 
        f | ((_namask & m_id ) >> 8) );
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

void AcdId::setRibbonNumber( UInt_t r )
{
    set_word( 0, m_id, r );
}

void AcdId::setRibbonOrientation( UInt_t orient) 
{
    set_word( 2, m_id, orient);
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
