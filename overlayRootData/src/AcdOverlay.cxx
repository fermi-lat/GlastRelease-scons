//									 
// The AcdOverlay class contains the information about a single ACD tile.
// This includes the PHA value, and above thresh information.  
//									 

#include "overlayRootData/AcdOverlay.h"
#include "Riostream.h"

ClassImp(AcdOverlay)

AcdOverlay::AcdOverlay() :
   m_volumeId(), m_acdId(), m_energyDep(0.), m_center(0.,0.,0.), m_status(0)
{
    return;
}

AcdOverlay::AcdOverlay(const VolumeIdentifier& volId) : 
            m_volumeId(volId),
            m_acdId(),
            m_energyDep(0.),
            m_center(0.,0.,0.),
            m_status(0)
{
    return;
}

AcdOverlay::AcdOverlay(const VolumeIdentifier& volId, 
                       const AcdId&            acdId, 
                       Double_t                energy, 
                       const TVector3&         position) : 
            m_volumeId(volId),
            m_acdId(acdId),
            m_energyDep(energy),
            m_center(position),
            m_status(0)
{
    return;
}

void AcdOverlay::Clear(Option_t *option) 
{
    TObject::Clear(option);
    return;
}

void AcdOverlay::Print(Option_t *option) const 
{
    using namespace std;
    TObject::Print(option);
    cout.precision(2);
    cout << "VolumeId: " << endl;
    m_volumeId.Print();
    cout << "AcdId: " << endl;
    m_acdId.Print();
    cout << "Energy = " << m_energyDep << endl;
    cout << "Position " << endl;
    m_center.Print();

    return;
}

Int_t AcdOverlay::Compare(const TObject *obj) const 
{
    const short base = 2;
    if (this == obj) return 0;
    if (AcdOverlay::Class() != obj->IsA()) return -1;
    UInt_t id_this = m_acdId.getId(base);
    UInt_t id_tile = ((AcdOverlay*)obj)->getAcdId().getId(base);
    if (id_this == id_tile)
	return 0;
    else
	return (id_this > id_tile) ? 1 : -1; 
}

Bool_t AcdOverlay::IsSortable() const {
  return kTRUE;
}

