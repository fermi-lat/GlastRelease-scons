#ifndef RootConvert_AcdReconConvert_H
#define RootConvert_AcdReconConvert_H 1

#include <Event/Recon/AcdRecon/AcdRecon.h>   
#include <reconRootData/AcdRecon.h>

#include <Event/Recon/AcdRecon/AcdReconV2.h>   
#include <reconRootData/AcdReconV2.h>

#include <Event/Recon/AcdRecon/AcdTkrPoint.h>   
#include <reconRootData/AcdTkrPoint.h>
#include <reconRootData/AcdTkrPointV2.h>

#include <Event/Recon/AcdRecon/AcdTkrHitPoca.h>
#include <reconRootData/AcdTkrHitPoca.h>
#include <reconRootData/AcdTkrHitPocaV2.h>

#include <Event/Recon/AcdRecon/AcdTkrGapPoca.h>
#include <reconRootData/AcdTkrGapPoca.h>
#include <reconRootData/AcdTkrGapPocaV2.h>

#include <Event/Recon/AcdRecon/AcdTkrIntersection.h>
#include <reconRootData/AcdTkrIntersection.h>

#include <Event/Recon/AcdRecon/AcdAssoc.h>
#include <reconRootData/AcdAssoc.h>

#include <Event/Recon/AcdRecon/AcdSplashVars.h>   
#include <reconRootData/AcdSplashVars.h>

#include "CLHEP/Matrix/SymMatrix.h"
#include <TMatrixDSym.h>


namespace RootPersistence {

  void convert( const Event::AcdRecon &, AcdRecon & ) ;
  void convert( const AcdRecon &, Event::AcdRecon & ) ; 

  void convert( const Event::AcdRecon &, AcdRecon & ) ;
  void convert( const AcdRecon &, Event::AcdRecon & ) ;

  void convert( const Event::AcdReconV2 &, AcdReconV2 & ) ;
  void convert( const AcdReconV2 &, Event::AcdReconV2 & ) ; 
  
  void convert( const Event::AcdHit&, AcdHit&);
  void convert( const AcdHit&, Event::AcdHit&);

  void convert( const Event::AcdEventTopology&, AcdEventTopology&);
  void convert( const AcdEventTopology&, Event::AcdEventTopology&);

  void convert( const Event::AcdAssoc&, AcdAssoc&);
  void convert( const AcdAssoc&, Event::AcdAssoc&);
  
  void convert( const Event::AcdTkrIntersection&, AcdTkrIntersection&);
  void convert( const AcdTkrIntersection&, Event::AcdTkrIntersection&);

  void convert( const Event::AcdTkrPoca&, AcdTkrPoca&);
  void convert( const AcdTkrPoca&, Event::AcdTkrPoca&);

  void convert( const Event::AcdTkrPoint&, AcdTkrPoint&);
  void convert( const AcdTkrPoint&, Event::AcdTkrPoint&);

  void convert( const Event::AcdTkrPoint&, AcdTkrPointV2&);
  void convert( const AcdTkrPointV2&, Event::AcdTkrPoint&);

  void convert( const Event::AcdTkrHitPoca&, AcdTkrHitPoca&);
  void convert( const AcdTkrHitPoca&, Event::AcdTkrHitPoca&);

  void convert( const Event::AcdTkrHitPoca&, AcdTkrHitPocaV2&);
  void convert( const AcdTkrHitPocaV2&, Event::AcdTkrHitPoca&);

  void convert( const Event::AcdTkrGapPoca&, AcdTkrGapPoca&);
  void convert( const AcdTkrGapPoca&, Event::AcdTkrGapPoca&);

  void convert( const Event::AcdTkrGapPoca&, AcdTkrGapPocaV2&);
  void convert( const AcdTkrGapPocaV2&, Event::AcdTkrGapPoca&);

  void convert( const Event::AcdSplashVars&, AcdSplashVars&);
  void convert( const AcdSplashVars&, Event::AcdSplashVars&);

  void convert( const Event::AcdTkrLocalCoords&, AcdTkrLocalCoordsV2&);
  void convert( const AcdTkrLocalCoordsV2&, Event::AcdTkrLocalCoords&);  

  void convert( const Event::AcdPocaData&, AcdPocaDataV2&);
  void convert( const AcdPocaDataV2&, Event::AcdPocaData&);  

  void convert( const CLHEP::HepSymMatrix&, TMatrixDSym& );
  void convert( const TMatrixDSym&, CLHEP::HepSymMatrix& );

  void convert( const CLHEP::HepSymMatrix&, TMatrixFSym& );
  void convert( const TMatrixFSym&, CLHEP::HepSymMatrix& );
}

#endif












