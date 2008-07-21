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

 @section gleam Gleam options
 For Gleam runs, append the following to basicOptions
 @verbatim
 // ------------------------
//    interleave background
//
ApplicationMgr.DLLs+= { "Interleave" };
Triggered.BranchMembers = {"TriggerAlg", "InterleaveAlg"}; //pseudo-sources do not trigger

FluxSvc.source_lib += {"$(INTERLEAVEXMLPATH)/source_library.xml"};
FluxAlg.sources += {"sampled_background"};
FluxAlg.sources += {"sampled_EarthAlbedo"};
InterleaveAlg.FileName= "$(INTERLEAVEXMLPATH)/exampleBackgroundSet.xml";
@endverbatim

  @section job JobOptions

@param Filename     set file name to pass to the backgroundSelection object  
@param TreeName     ["MeritTuple"] name of tree
@param DisableList List of branches to disable from being copied. Default, which can be
appended to or replaced, is {"EvtElapsedTime","EvtLiveTime","Pt*","FT1*","CT*"}. These are the 
quantities that will *not* be copied from the background tuple, and must be set according to the
current run.
@param MapName     ["interleave_map"] If set, a ROOT TTree with this name will be added to the
output file, with a correspondence of run/event numbers from the sampled file to the output

  <hr>
  \section notes release notes
  release.notes
  \section requirements requirements
  \include requirements

*/

