/** @mainpage package mcRootData

  This package defines the ROOT classes which mirror the GlastEvent MonteCarlo
  TDS classes.  The primary ROOT Monte Carlo class is McEvent.  
  
  @image html mcTree.gif
    
  Each McEvent contains:

  - Event Id
  - Run Number
  - Collection of McParticles
    -# particle id
    -# status flags
    -# initial position
    -# final position
    -# initial momentum
    -# final momentum
    -# reference (TRef) to mother McParticle
    -# collection of McParticle daughters
  - Collection of McPositionHit
    -# monte carlo particle id
    -# volume id
    -# reference (TRef) to the origin McParticle
    -# entry position
    -# exit position
    -# deposited energy
    -# particle energy
    -# time of flight
    -# direction cosine - calculated at run time
  - Collection of McIntgratingHit
    -# volume id
    -# reference (TRef) to an McParticle
    -# status flags
    -# first and second moments
    -# map of McParticle and energies

  @attention In ROOT 3.02.07, it is not possible to access TRefs stored in a
  ROOT file from compiled code.  This bug has been fixed in the development 
  version of ROOT.  For now, one must use ROOT interactively to gain access to 
  TRefs stored in a ROOT file.

  <hr>
  @section notes release notes
  release.notes
  <hr>
  @section requirements requirements
  @verbinclude requirements

*/

