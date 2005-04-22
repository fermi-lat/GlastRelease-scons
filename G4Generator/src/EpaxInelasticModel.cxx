#include "EpaxInelasticModel.hh"
#include "Randomize.hh"

G4HadFinalState *EpaxInelasticModel::ApplyYourself
                  ( const G4HadProjectile &aTrack , G4Nucleus &targetNucleus )
{ 
//   G4cout << "In Epax!" << G4endl;
   IsApplicable ( aTrack );
   G4float fSigmaOfReaction;
   G4int ProjA = (G4int) ( aTrack.GetDefinition()->GetPDGMass() / amu_c2 + 0.5 );
   G4int ProjZ = (G4int) ( aTrack.GetDefinition()->GetPDGCharge() / eplus + 0.5 );
   //G4cout << "ProjZ:" << ProjZ << " ProjA:" << ProjA << G4endl;
   G4int TargA = (G4int) ( targetNucleus.GetN() + 0.5 );
   G4int TargZ = (G4int) ( targetNucleus.GetZ() + 0.5 );
   G4double barrier;
   const G4LorentzVector proj4Momentum = aTrack.Get4Momentum();
   G4ThreeVector incidentVector_ini = proj4Momentum.getV();
   G4double beta=proj4Momentum.beta();
   G4ThreeVector incidentVector = incidentVector_ini.unit();
   G4ThreeVector p00=(1./ProjA) * incidentVector_ini; 
   const G4Material *TargetMaterial = aTrack.GetMaterial();
   G4int NumberOfElement = TargetMaterial->GetNumberOfElements();
   const G4Element *TargetElement = 0;
   G4double ZofElement;
   G4double NofElement; 
   G4int i; 
   for ( i = 0 ; i < NumberOfElement ; i ++ ) 
   {
      ZofElement = TargetMaterial->GetElement(i)->GetZ();
      NofElement = TargetMaterial->GetElement(i)->GetN();
      if ( abs ( ZofElement - targetNucleus.GetZ() ) < 1.0 )
         TargetElement = TargetMaterial->GetElement(i);
   }
   G4double Sigma;
   G4ParticleDefinition* pd = (G4ParticleDefinition*)aTrack.GetDefinition();
   const G4DynamicParticle* originalIncident = new G4DynamicParticle (pd, proj4Momentum);   
   if ( (aTrack.GetKineticEnergy()/ProjA)   < 1000 )
     Sigma = TCS.GetCrossSection(originalIncident, TargetElement,0);
   else
     Sigma = SCS.GetCrossSection(originalIncident, TargetElement,0);
   delete originalIncident;
   fSigmaOfReaction = ( G4float ) ( Sigma / barn );
   theParticleChange.SetEnergyChange( 0.0 );
   theParticleChange.SetStatusChange(stopAndKill);
   G4DynamicParticle *nDP; 
   nDP = new G4DynamicParticle();  
   G4PT = G4ParticleTable::GetParticleTable();
   G4ThreeVector SecondaryMomentum;
   G4int mult,Z,A;
   //  
   select_ion(fSigmaOfReaction,TargA,TargZ,ProjA,ProjZ,&A,&Z); 
   //G4cout << TargA << " " << TargZ  << " " << Z << " " << A << G4endl;
   if (Z != 1 /*&& Z !=4*/) {
      G4double ExcitationEnergy = ( G4double ) 0 * MeV ;
        nDP -> SetDefinition ( G4PT -> GetIon ( Z , A , ExcitationEnergy ) );
        SecondaryMomentum=A*p00;
	nDP -> SetMomentum ( SecondaryMomentum );
         theParticleChange.AddSecondary( nDP );
	 mult= (ProjA-A)/2; }
   else { 
     // desintegration
       if (Z==1) mult= (ProjA)/2;
     // breakup of 8Be into 2 alphas 
       /*  else {
           nDP -> SetDefinition (  G4Alpha::Alpha() );
           SecondaryMomentum=A*p00;
 	   nDP -> SetMomentum ( SecondaryMomentum );
           theParticleChange.AddSecondary( nDP );
	   nDP = new G4DynamicParticle();
           nDP -> SetDefinition (  G4Alpha::Alpha() );
	   nDP -> SetMomentum (SecondaryMomentum);
           theParticleChange.AddSecondary( nDP );
	   mult=(ProjA-4)/2; }*/
   }
   G4double fact_p[3]={1.,.7,0.};
   if (TargA < 30) fact_p[0]=.3;
   G4double fact_mult[3]={1.,1.,2.};
   G4double fact_mult_a[3]={.2,.2,0.};
   G4double fact_T[3]={3.,5.,3.};
   G4double A_s[3], Z_s[3];
   G4int A_p[3]={1,1,4};   
   A_s[0]=ProjA;
   Z_s[0]=ProjZ;
   A_s[1]=2*(ProjA-A);
   Z_s[1]=2*(ProjZ-Z);
   A_s[2]=TargA;
   Z_s[2]=TargZ;
   for (G4int ipart=0; ipart<3; ipart++) { // loop over neutrons,protons & alphas
          for (G4int icont=0; icont<3; icont++) { 
	     barrier=0;
             if (!ipart) {
		G4double r0=1.3*(pow(A_s[icont],0.333)+A_p[ipart]); 
	        barrier=1.44*Z_s[icont]*ipart/r0; } // Z_p = ipart 
             G4int mult_s=(G4int) mult*fact_mult[icont];
             if (ipart==2) mult_s *= G4UniformRand()*fact_mult_a[icont]; 
             G4double T=sqrt(mult_s*200./ProjA)*fact_T[icont];  
             for (int i = 0 ; i < mult_s ; i ++ ) {
               nDP = new G4DynamicParticle();  
               if (ipart==0) nDP -> SetDefinition ( G4Neutron::Neutron() );  
               if (ipart==1) nDP -> SetDefinition ( G4Proton::Proton() );
               if (ipart==2) nDP -> SetDefinition ( G4Alpha::Alpha() ); 
               G4double M_p= nDP->GetMass();
		// isotropic distribution
                 G4double cost=(1-2*G4UniformRand());
                 G4double sint=sqrt(1-cost*cost);
                 G4double phi=2*pi*G4UniformRand();
                 // maxwellian spectrum
                 G4double chi1= G4UniformRand(); 
                 G4double chi2= G4UniformRand(); 
                 G4double chi3= G4UniformRand() ;
                 G4double ep=T *(-log(chi1)-log(chi2)*cos(pi * chi3 /2)*cos(pi * chi3 /2) );    //MeV
                  G4double E_tot=ep+barrier+M_p;
		  G4double pp=sqrt(E_tot*E_tot-M_p*M_p); // MeV
		  // G4cout << icont <<" " << ipart << " " << pp  << " " << T << G4endl;        
		  G4LorentzVector pcm(pp*cost,pp*sint*cos(phi),pp*sint*sin(phi),E_tot); 
		 //Lorentz transform to lab frame
                 SecondaryMomentum=(pcm.boost(incidentVector,fact_p[icont]*beta)).getV();		 
                 nDP -> SetMomentum ( SecondaryMomentum );
                 theParticleChange.AddSecondary( nDP );
	     }
	  }
       }
   return &theParticleChange;
}

