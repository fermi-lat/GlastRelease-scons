/** @mainpage package mcRootData

  This package defines the ROOT classes which mirror the GlastEvent MonteCarlo
  TDS classes.  The primary ROOT Monte Carlo class is McEvent.  Each McEvent 
  contains:

  - Event Id
  - Run Number
  - Collection of McParticles
    -# particle id
    -# status flags
    -# initial position - calculated at run time
    -# final position
    -# initial momentum
    -# final momentum
    -# reference (TRef) to mother McParticle
    -# collection of daughter McParticles
  - Collection of McPositionHits
    -# volume id
    -# reference (TRef) to the McParticle
    -# reference (TRef) to the origin McParticle
    -# entry position
    -# exit position
    -# deposited energy
    -# particle energy
    -# time of flight
    -# direction cosine - calculated at run time
  - Collection of McIntgratinghits
    -# volume id
    -# reference (TRef) to an McParticle
    -# status flags
    -# first and second moments
    -# map of McParticle and energies

  @attention In ROOT 3.02.07, it is not possible to access TRefs stored in a
  ROOT file from compiled code.  This bug has been fixed in the development 
  version of ROOT.  For now, one must use ROOT interactively to gain access to 
  TRefs stored to a ROOT file.

  <hr>
  @section notes release notes
  release.notes
  @section requirements requirements
  @include requirements

*/

