/** @mainpage package Gleam
@author Toby Burnett

@section intro Introduction
This package defines the Gleam application

  @section used Used Packages for building

The requirements file lists explicit versions for all packages required to build 
the application's component dll's. See the requirements file. 

  @section build Build Targets

It also builds (redundantly) two basic executables from the source files in GlastPolicy and GuiSvc.

  - Gleam Gui-based. Run from the requirements files in the src folder
  - test_Gleam Non gui. Run from src/test/jobOptions.txt


@section joboptions job options files

The following files are in the <a href="../../src">src folder</a>:

 <ul>
   <li> <a href ="../../src/basicOptions.txt">basicOptions.txt</a>

  Meant to be #included by all other options files. 
  All basic services are defined. Defines the following Sequencers:
   <ul>
     <li> Top - top level: this is the only Algorithm  in the ApplicationMgr.TopAlg list. 
     It is initialized to other Sequencers :Generator, Digitization, EventDisplay, Triggered.
     <li> Generator 
     <li> Digitization
     <li> EventDisplay
     <li> Triggered 
     <li> Reconstruction
     <li> RecoDisplay
     <li> Output
   </ul>
   The following DLL's are defined:
   <ul>
     <li> GaudiAlg
     <li> GlastSvc
     <li> G4Generator
     <li> TkrDigi
     <li> CalDigi
     <li> AcdDigi
     <li> Trigger
     <li> Recon
     <li> CalRecon
     <li> TkrRecon
     <li> AcdRecon
     <li> RootIo
   </ul>
  <li> <a href ="../../src/guiOptions.txt">guiOptions.txt</a>

  Meant to be included after basicOptions.  Appends GuiSvc to the DLLs, and to the ExtSvc list. 
  Sets ApplicationMgr.Runnable to GuiSvc to allow control of the event loop.

  <li> <a href ="../../src/default.txt">default.txt</a> 

  For standard options. No GUI, merit, default source.

  
  <li> <a href ="../../src/gamma_1_gev_normal.txt">gamma_1_gev_normal.txt</a>

  Not functional. TODO: fix, or remove.

  <li> <a href ="../../src/jobOptions.txt">jobOptions.txt</a>

  Standard job: sets up GUi for display, in pause (interactive) mode. Defines the muon pencil source.

  <li> <a href ="../../src/readdigi.txt">readdigi.txt</a>

  Really an example of how to modify the basicOptions to configue the application to read from a 
  root file containing digitization information. Assumes that a file "digi.root" exists in the default 
  directory. Also, job crashses at the end due to error returned by digiRootReaderAlg.

</ul>

@section systests System Tests
The folder <a href="../../src/systests">systests</a> Contains system tests.
<hr>
@section  jobOptions jobOptions
Since there is no code, there are no specific job options
$Heading:$
<hr>
@section notes release notes
@include release.notes
@section requirements requirements
@include requirements

*/
