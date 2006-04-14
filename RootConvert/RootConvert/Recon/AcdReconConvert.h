#ifndef RootConvert_AcdReconConvert_H
#define RootConvert_AcdReconConvert_H 1

#include <Event/Recon/AcdRecon/AcdRecon.h>   
#include <reconRootData/AcdRecon.h>

#include <Event/Recon/AcdRecon/AcdTkrPoint.h>   
#include <reconRootData/AcdTkrPoint.h>

#include <Event/Recon/AcdRecon/AcdTkrHitPoca.h>
#include <reconRootData/AcdTkrHitPoca.h>

#include <Event/Recon/AcdRecon/AcdTkrGapPoca.h>
#include <reconRootData/AcdTkrGapPoca.h>

namespace RootPersistence {

  void convert( const Event::AcdRecon &, AcdRecon & ) ;
  void convert( const AcdRecon &, Event::AcdRecon & ) ; 

  void convert( const Event::AcdRecon &, AcdRecon & ) ;
  void convert( const AcdRecon &, Event::AcdRecon & ) ;
  
  void convert( const Event::AcdHit&, AcdHit&);
  void convert( const AcdHit&, Event::AcdHit&);

  void convert( const Event::AcdTkrIntersection&, AcdTkrIntersection&);
  void convert( const AcdTkrIntersection&, Event::AcdTkrIntersection&);

  void convert( const Event::AcdTkrPoca&, AcdTkrPoca&);
  void convert( const AcdTkrPoca&, Event::AcdTkrPoca&);

  void convert( const Event::AcdTkrPoint&, AcdTkrPoint&);
  void convert( const AcdTkrPoint&, Event::AcdTkrPoint&);

  void convert( const Event::AcdTkrHitPoca&, AcdTkrHitPoca&);
  void convert( const AcdTkrHitPoca&, Event::AcdTkrHitPoca&);

  void convert( const Event::AcdTkrGapPoca&, AcdTkrGapPoca&);
  void convert( const AcdTkrGapPoca&, Event::AcdTkrGapPoca&);

}

#endif












