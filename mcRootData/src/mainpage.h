/** @mainpage package mcRootData

  This package defines the ROOT classes which mirror the GlastEvent MonteCarlo
  TDS classes.  The primary ROOT Monte Carlo class is McEvent.  Each McEvent 
  contains:

  - Event Id
  - Run Number
  - List of McParticles
    -# Reference (TRef) to mother McParticle
  - List of McPositionHits
    -# Reference (TRef) to the McParticle
  - List of McIntgratinghits
    -# Reference (TRef) to an McParticle

  <hr>
  @section notes release notes
  release.notes
  @section requirements requirements
  @include requirements

*/

