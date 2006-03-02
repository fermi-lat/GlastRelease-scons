//                                                                       
// The L1T class contains the information about the Level 1 Trigger count
// for an event.  
//                                                                       

#include "digiRootData/Tem.h"
#include "Riostream.h"

ClassImp(Tem)


Tem::Tem() 
{
    Clear();
}

Tem::Tem(UShort_t t, const ErrorData& err) {
    Clear();
    init(t, err);
}

Tem::Tem(const Tem& t):TObject(t) {
    Clear();
    m_tower = t.m_tower;
    m_error = t.m_error;
}

void Tem::init(UShort_t t, const ErrorData& err) {
    m_tower = t;
    m_error = err;
}


void Tem::Clear(Option_t *option) {
    m_tower = 0;
    m_error.Clear(option);
}

void Tem::Print(Option_t *option) const {
    using namespace std;
    TObject::Print(option);
    cout.precision(2);
    cout << "Tower: " << m_tower << endl;
    m_error.Print();
   
}


