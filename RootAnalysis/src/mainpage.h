
/** @mainpage package RootAnalysis

  @section intro Introduction

  This package contains all GLAST ROOT macros that are supported by the GLAST 
  ground software group.  The routines provided, are meant to provide some 
  basic, general, functions that most GLAST ROOT users may find useful.  
  To process GLAST ROOT data, it is not necessary to use this package, but it 
  may make it easier.

  @section setup Setup the ROOT environment for GLAST

  RootAnalysis now comes with a simple .rootrc.  This is the file used by ROOT
  to setup your ROOT environment upon startup.  
  First, creat a ROOTANALYSIS environment variable that is set to the directory 
  path, where you have the RootAnalysis package installed. 

  On UNIX:  setenv ROOTANALYSIS $HOME/RootAnalysis/vnrm

  On Windows: 
  - In Control Panel, Click on System
  - Click on the Environment tab 
  - Define a new variable called ROOTANALYSIS
  Set its value equal to the RootAnalysis directory on your system, 
  i.e. D:/glast/RootAnalysis/v2r2 

  To use the RootAnalysis package, the .rootrc file must be in the directory 
  from where you launch ROOT.  Next do one of the following:

  Start ROOT from within the RootAnalysis directory
  OR
  Copy the .rootrc file to your local area, where you will start ROOT. 

  @section cmt Using RootAnalysis as a regular CMT package
  It is also possible to download RootAnalysis directly from the CVS 
  repository.  In this case, you must make sure to also install the 
  mcRootData, digiRootData, and reconRootData packages on your system.
  In this case, the RootAnalysis environment variable will be set up for you.
  Follow the typical directions for setup using CMT.

  Note that RootAnalysis now contains a test routine that runs RootTreeAnalysis
  from a standalone main program.  By default, RootTreeAnalysis reads in 3
  ROOT files that come with the standard RootAnalysis distribution.  This test
  routine, can be run on other ROOT files, by providing input parameters:
  test_RootAnalysis.exe myMcRootFile.root myDigiRootFile.root myReconRootFile.root

  @section loop Event Loop Processing

  \b RootTreeAnalysis
  
    Our generic event loop macro.  Handles digi, reconstruction, and/or 
    monte carlo ROOT files.  Output is a new ROOT file which contains 
    user defined histograms.  Users can provide either file names, or 
    TChains of files for multi-file processing.

  \b NtupleAnalysis
  
    This class handles any ntuple file, and produces an output PostScript 
    file containing plots of any and all entries in the ntuple - and can 
    apply an optional cut to the data as well before creating the plots.
    One can pre-set the number of columns and rows of the PS pages as well - 
    to determine how many plots to display per page

  @section low Low Level Utilities
  Some low level utilities are made available in the RootAnalysis/utilityMacro 
  directory.  Each macro file contains some description of the macro's 
  function and its usage.

  \b chainTrees.cxx
  
    This small macro accepts a list of files names, the path and name of the 
    TTree and returns a TChain.  All files must contain TTrees of the same 
    type i.e. ntuple, digi, recon, or Monte Carlo.  A TChain is useful if you 
    want to analyze multiple files at once.  The TChain created can then be 
    passed to RootTreeAnalysis for event loop processing.

  \b chainAll.cxx
  
    Similar to chainTrees, this macro chains all Trees located in ROOT files in 
    a specified directory.  ROOT files are identified by having the *.root 
    extension - this macro will ignore all other files.

  \b copyTree.cxx
  
    Copies a TTree from one file into a new file, allowing for a new directory to be 
    created to store the TTree in the new file.

  \b createEventList.cxx
  
    This macro provides a generic mechanism to create a ROOT TEventList.  
    Event lists are useful, in that they restrict the entries in a TTree that 
    will be processed.   A TEventList is created, by applying a cut on the data.
    There are 2 versions of this macro contained in the same file, one for a 
    regular TFile and the other for a TChain (a chain of files):

  \b pruneTree.cxx
  
    This macro creates a copy of a TTree containing a truncated eventlist 
    based upon user supplied cuts.

  <hr>
  @section notes release notes
  release.notes
  @section requirements requirements
  @include requirements

*/

