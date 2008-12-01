#include "overlayRootData/CalOverlay.h"
#include "Riostream.h"

ClassImp(CalOverlay)

//
// Class:  CalOverlay        
// Description: Digitizations for Cal
// Based on the CalOverlay TDS class in GlastEvent written by J. Eric Grove.
//              
// Author:  Heather Kelly
//

CalOverlay::CalOverlay() : 
    m_position(0.,0.,0.), m_energy(0.), m_status(0)
{
    return;
}

CalOverlay::~CalOverlay() { 
    return;
}

void CalOverlay::initialize(const CalXtalId &id, const TVector3& position, const Double_t& energy, const UInt_t& status) 
{
    m_xtalId = id;
    m_position = position;
    m_energy   = energy;
    m_status   = status;
}

void CalOverlay::Clear(Option_t *option) {
    m_xtalId.Clear();
    m_position = TVector3(0.,0.,0.);
    m_energy   = 0.;
    m_status   = 0;
}

void CalOverlay::Print(Option_t *option) const 
{
    TObject::Print(option);
    cout.precision(2);
    std::cout << "Id: " << endl;
    m_xtalId.Print();

    std::cout << "Energy: " << m_energy << ", position: ";
    m_position.Print();
    std::cout << std::endl;
}
