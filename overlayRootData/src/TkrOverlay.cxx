#include "overlayRootData/TkrOverlay.h"
#include "Riostream.h"

//#include "DigiObjectManager.h"

ClassImp(TkrOverlay)

TkrOverlay::TkrOverlay() {
    Clear();
}

TkrOverlay::~TkrOverlay (){
    Clear();
}

void TkrOverlay::Clear(Option_t *option) {
    m_bilayer = 0;
    m_tower = 0;
    m_tot[0] = 0;
    m_tot[1] = 0;
    m_lastController0Strip = -1;
    m_hitCol.clear();
}

void TkrOverlay::Print(Option_t *option) const {
    using namespace std;
    TObject::Print(option);
    cout.precision(2);
    cout << "BiLayer: " << m_bilayer << endl;
    m_tower.Print();
    cout << "ToT: " << m_tot[0] << " " << m_tot[1]
        << "  LastController0Strip: " << m_lastController0Strip << endl;
    cout << "NumHits: " << getNumHits() << endl;
}

void TkrOverlay::initialize(Int_t l, GlastAxis::axis v, TowerId t, Int_t* tot)
{
    m_bilayer = l;
    m_tower = t;
    m_view = v;
    m_tot[0] = tot[0];
    m_tot[1] = tot[1];
    
    m_lastController0Strip =  -1;    
}

Int_t TkrOverlay::getHit(UInt_t i) const
{
    return (i < m_hitCol.size() ? m_hitCol[i] : -1);
}

Int_t TkrOverlay::getStrip(UInt_t i) const
{
    return (i < m_hitCol.size() ? m_hitCol[i] : -1);
}

Int_t TkrOverlay::getToT(UInt_t i) const {
    if (i < 2) return m_tot[i];
    return -1;
}

//! Get ToT for a given strip
Int_t TkrOverlay::getToTForStrip(UInt_t strip) const
{ 
    // Purpose and Method:  Returns the ToT for a given strip by
    //  determining which controller this strip is associated with.
    //  We then return the ToT corresponding to that controller.
    if (m_lastController0Strip == -1) return -1;
    return m_tot[(Int_t(strip)<=m_lastController0Strip ? 0 : 1 )]; 
}


void TkrOverlay::addC0Hit( UInt_t strip ) 
{
    // Purpose and Method:  Add a controller 0 strip to the hit list
    //   and keeps track of highest strip associated with controller 0
    m_hitCol.push_back(strip);  
    if (m_lastController0Strip < (Short_t) strip) m_lastController0Strip = strip;  
}

//! Add a controller 1 hit to the hit list
void TkrOverlay::addC1Hit( UInt_t strip ) 
{
    m_hitCol.push_back(strip);
}


//void* TkrOverlay::operator new(size_t size)
//{
//    TkrOverlay* temp = DigiObjectManager::getPointer()->getNewTkrOverlay();
//
//    // Since we recycle, make sure these member functions are cleared
//    temp->m_hitCol.clear();
//
//    return temp;
//}

//void* TkrOverlay::operator new(size_t size, void* vp)
//{
//    return vp;
//}

//void TkrOverlay::operator delete(void* p)
//{
//    // Since we let DigiObjectManager handle memory, nothing to do here
//    return;
//}

