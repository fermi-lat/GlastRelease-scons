// Classe che gestisce gli hit sul layer
#ifndef LayerHit_h
#define LayerHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class LayerHit : public G4VHit
{
  public:

      LayerHit();
     ~LayerHit();
      LayerHit(const LayerHit&);
      const LayerHit& operator=(const LayerHit&);
      int operator==(const LayerHit&) const;

      inline void* operator new(size_t);
      inline void  operator delete(void*);

      void Draw();
      void Print();
      
  public:
  
      void AddSil(G4double de, G4double dl) {EdepSil += de; TrackLengthSil += dl;};
                 
      G4double GetEdepSil()     { return EdepSil; };
      G4double GetTrakSil()     { return TrackLengthSil; };
     
  private:
  
      G4double EdepSil, TrackLengthSil;
      
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

typedef G4THitsCollection<LayerHit> LayerHitsCollection;

extern G4Allocator<LayerHit> LayerHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void* LayerHit::operator new(size_t)
{
  void* aHit;
  aHit = (void*) LayerHitAllocator.MallocSingle();
  return aHit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void LayerHit::operator delete(void* aHit)
{
  LayerHitAllocator.FreeSingle((LayerHit*) aHit);
}

#endif


