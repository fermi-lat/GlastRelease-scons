#include "reconRootData/Recon.h"

ClassImp(Recon)

// default constructor does not create
// any objects
Recon::Recon() {
}

Recon::~Recon() {
    Clean();
}

void Recon::Clean() {
    
    m_recFlags.Clear();
    m_cal.Clean();
    m_tkr.Clean();
}

void Recon::Create() {

    m_tkr.Create();
    m_cal.Create();

}
