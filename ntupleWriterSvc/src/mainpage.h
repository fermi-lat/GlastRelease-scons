// $Header$
// Mainpage for doxygen

/** @mainpage package ntupleWriterSvc
 * @author Heather Kelly, T Burnett
 *
 * @section intro Introduction
 
 * Defines the service RootTupleSvc, which implements INTupleWriterSvc. It allows for multiple trees, and has only a setup phase:
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
 ...
    // each event
    m_rootTupleSvc->storeRowFlag(true); 

@endverbatim

 The properties are listed in the following table:

 @param filename [RootTupleSvc.root]  The ROOT filename	
 @param checksumfilename [""]
 @param treename	[1] Default tree name to use if the "tupleName" argument is empty
 @param title	["Glast tuple"] Title for all TTree's
 @param defaultStoreFlag [false] Store the row at the end of an event loop cycle. If false, client must execute storeRowFlag(true) to save all trres in  the current event. Or, the client can execute storeRowFlag(treename, true) to save the specific tree.	false
 @param AutoSave [100000]	 Number of bytes to save at a time (ROOT default is 10000000)

 * <hr>
 * @section notes release notes
 * release.notes
 * @section requirements requirements
 * @verbinclude requirements
*/

