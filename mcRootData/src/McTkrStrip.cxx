#include <mcRootData/McTkrStrip.h>
#include <commonRootData/RootDataUtil.h>
#include "Riostream.h"

#include "McObjectManager.h"

ClassImp(McTkrStrip)

McTkrStrip::McTkrStrip(const VolumeIdentifier& id, UInt_t strip, Double_t e,
         Bool_t noise, Double_t deltaX, Double_t deltaY, TRefArray* hitList)
        : m_planeId(id)
        , m_strip(strip)
        , m_energy(e)
        , m_noise(noise)
        , m_deltaX(deltaX)
        , m_deltaY(deltaY) 
{
    m_hits.Clear();
    if (hitList) 
    {
        for(int idx=0; idx < hitList->GetEntries(); idx++)
            m_hits.Add(hitList->At(idx));
    }
};

McTkrStrip::McTkrStrip()
        : m_strip(0)
        , m_energy(0)
        , m_noise(false)
        , m_deltaX(0.)
        , m_deltaY(0.) 
{
    m_hits.Clear();
};

McTkrStrip::~McTkrStrip()
{
    Clear();
}
    
void McTkrStrip::initialize(const VolumeIdentifier& id, UInt_t strip, Double_t e,
         Bool_t noise, Double_t deltaX, Double_t deltaY, TRefArray* hitList)
{
    m_planeId = id;
    m_strip   = strip;
    m_energy  = e;
    m_noise   = noise;
    m_deltaX  = deltaX;
    m_deltaY  = deltaY;

    m_hits.Clear();
    if (hitList) 
    {
        for(int idx=0; idx < hitList->GetEntries(); idx++)
            m_hits.Add(hitList->At(idx));
    }
};

void McTkrStrip::Clear(Option_t *)
{
    m_hits.Clear();
}

void McTkrStrip::Print(Option_t *option) const {
    using namespace std;
    TObject::Print(option);
    m_planeId.Print(option);
    std::cout.precision(2);
    std::cout << "Strip #" << m_strip << ", deposited energy=" << m_energy << std::endl;
}
