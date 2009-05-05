// $Header$
// Mainpage for doxygen

/** @mainpage package ntupleWriterSvc
 * @author Heather Kelly, T Burnett
 *
 * @section intro Introduction
 
 * Defines the service RootTupleSvc, which implements INTupleWriterSvc. It allows for multiple trees and multiple files, and has only a setup phase:
 one defines the tuple at initialization time, with pointers to the tuple values. 

 <br> Usage example: (see writeJunkAlg for a complete example)

 @verbatim
 #include "ntupleWriterSvc/INTupleWriterSvc.h"
 ...
    // declarations
    INTupleWriterSvc *m_rootTupleSvc;
    double test;
 ...
    // during setup
    sc = service("RootTupleSvc", m_rootTupleSvc);
    m_rootTupleSvc->addItem("","test", &test);
    or
    m_rootTupleSvc->addItem("","test", &test, "myFile.root");
 ...
    // each event
    m_rootTupleSvc->storeRowFlag(true); 

@endverbatim

 <hr>
 * @section jobOptions jobOptions
 * @param RootTupleSvc.filename 
 * Default RootTupleSvc.root
 * Name of the output ROOT file that will contain the ROOT ntuple
 * @param RootTupleSvc.treename
 * Default "1"
 * Logical name of the ROOT TTree that will contain the ntuple data
 * @param RootTupleSvc.title
 * Default "Glast tuple"
 * Title of the ROOT TTree
 * @param RootTupleSvc.defaultStoreFlag
 * Default false
 * Store the row at the end of an event loop cycle. If false, client must 
 * execute storeRowFlag(true) to save all trres in  the current event. Or, the 
 * client can execute storeRowFlag(treename, true) to save the specific tree.
 * @param RootTupleSvc.AutoSave 
 * Default 100000
 * In Bytes, denoting the size the ntuple must reach before triggering a true
 * write to disk.  ROOT's default is 10000000
 * @param RootTupleSvc.RejectIfBad
 * Default true
 * if set, tuple entries containing any non-finite values are not written
 * @param RootTupleSvc.StartingIndex
 * Default 0
 * Used for input ROOT tuples to denote starting index to read
 * @param RootTupleSvc.JobInfoTreeName 
 * Default ["jobinfo"]
 * Name of the tree to contain the jobinfo data
 * @param RootTupleSvc.MeritVersion
 * Default 0 which results in using the initialized value set in RootTupleSvc
 * @param RootTupleSvc.JobInfo  
 * Default "" (empty string)
 * string of the form "a=99, b=100" of names, values to add to the jobinfo tree
 * some use this string to store the tag of GR, such as 
 * "version=13, revision=2, patch=3"
 * @param RootTupleSvc.IncludeBranchList
 * Default "" (empty string)
 * A list of merit branch names to read when reading an input set of merit files
 * this list may use wildcards, and must conform to the case-sensitive names
 * of the actual branches in the ROOT file
 * @param RootTupleSvc.ExcludeBranchList
 * Default "" (empty string)
 * A list of merit branch names to exclude when reading an input set of merit 
 * files this list may use wildcards, and must conform to the case-sensitive 
 * names of the actual branches in the ROOT file
 * <hr>
 * @section notes release notes
 * release.notes
 * @section requirements requirements
 * @verbinclude requirements
*/

