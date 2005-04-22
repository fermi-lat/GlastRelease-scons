//
// JQMD Inelastic Model Geant4-JQMD
// T. Koi NIRS iSRL 01-Jul-2002
// T. Koi SLAC SCS  tkoi@slac.stanford.edu
//
// Class Description
// Final state production model for GenericIon  inelastic scattering 
// below 3.5 GeV/Nucleon 
// To be used in your physics list in case you need this physics.
// In this case you want to register an object of this class with 
// the corresponding process.
// Class Description 
//
// 15-Oct-2003 This Version written by  T. Koi
// 15-Dec-2003 Adapete to Geant4 Ver6. by T. Koi
//             Change IsApplicapble
//
 
#ifndef EpaxInelasticModel_h
#define EpaxInelasticModel_h 1

#include "G4InelasticInteraction.hh"
#include "G4EpaxFragmentCrossSection.hh"

#include "TripathiCrossSection.hh"
#include "G4IonsSihverCrossSection.hh"

#include "G4ParticleTable.hh"


class EpaxInelasticModel : public G4InelasticInteraction
{
   public:
    
      EpaxInelasticModel() : G4InelasticInteraction()
      {
         SetMinEnergy( 0.0 ); // This value is temporary.
                         // Strongly suggested change several tens MeV
         SetMaxEnergy ( 100 * 240 * GeV ); // 100 GeV/nucleon
      }
    
      ~EpaxInelasticModel() { }
    
//      G4VParticleChange *ApplyYourself( const G4Track &aTrack,
      G4HadFinalState *ApplyYourself( const G4HadProjectile &aTrack,
                                      G4Nucleus &targetNucleus );

//    Check Energy of Projetile
      void IsApplicable ( const G4HadProjectile &aTrack )
      {
        /* G4int ProjA = (G4int) ( aTrack.GetDefinition()->GetPDGMass() / amu_c2 + 0.5 );  
         //const G4DynamicParticle *DP = aTrack.GetDynamicParticle();
         //G4double ProjKe = DP->GetKineticEnergy();
         G4double ProjKe =  aTrack.GetKineticEnergy();
	*/
      }
  void select_ion(G4double sigma_tot, G4int At, G4int Zt, G4int Ap, G4int Zp, G4int* A0, G4int* Z0){
 
G4EpaxFragmentCrossSection i;
G4double st,sigr[30],res;
 G4int A,Z;
 sigr[3]=0;
   st=0;
   for (G4int izf=4; izf<Zp; izf++){
     sigr[izf]=0;
     res=0;
     G4int amin=2*izf-7;
     if (amin < izf+2) amin=izf+2;
          for (G4int ia=amin; ia<2*izf+8; ia++){
       if (ia > Ap || (ia-izf) > (Ap-Zp)) break; 
        A=ia;
        Z=izf; 
       res+=i.doit(Ap,Zp,At,Zt,A,Z);       
}
     st+=res;
     sigr[izf]=sigr[izf-1]+res;
} 
 G4double xran=G4UniformRand();
 if (Zp > 4 && xran < st/sigma_tot)  {xran=G4UniformRand();
              for (Z=4;Z<Zp;Z++){ if (xran < sigr[Z]/st) break;}
              *A0=2*Z;
              *Z0=Z;
 } else {*Z0=1;*A0=1; } 
}
   private:
      G4ParticleTable* G4PT;
      TripathiCrossSection TCS;
      G4IonsSihverCrossSection SCS;

      G4RotationMatrix RotMat;
      G4RotationMatrix inverseRotMat;
      G4ThreeVector RotAxis;
      G4double RotAngle;
};
 
#endif
 
