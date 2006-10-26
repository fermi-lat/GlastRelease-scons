// $Header$
// Mainpage for doxygen

/*! \mainpage package Interleave

   \author Toby Burnett <tburnett@u.washington.edu>
   

  This package manages the insertion of interleaved background events into a Gleam simulation.

  - InterleaveAlg   Gaudi algorithm that manages the insertion
  - SampledBackground Spectrum to insert flag particles at the trigger rate
  - BackgroundSelection Manage selection of background events from external database

  It is expected to be invoked by a branch from the Trigger Sequencer, as in the following job options invoked by the test program, see test.cxx


  @verbatim
  //the top sequence loop
Top.Members={"Sequencer/Generator"};

Generator.Members     = {"FluxAlg", "PointInfoAlg", "Sequencer/Trigger"};

Trigger.Members       = {"TriggerAlg"};
Trigger.BranchMembers = {"TriggerAlg", "TestAlg", "InterleaveAlg"};
 @endverbatim

  @section job JobOptions

@param RootFile     set root file name to pass to the backgroundSelection object  
@param TreeName     ["MeritTuple"] name of tree

  <hr>
  \section notes release notes
  release.notes
  \section requirements requirements
  \include requirements

*/

