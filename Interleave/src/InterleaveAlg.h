/** @file InterleaveAlg.h

@brief declaration  of the class InterleaveAlg

$Header$

*/
#ifndef InterleaveAlg_h
#define InterleaveAlg_h

#include "BackgroundSelection.h"

#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/Property.h"
class BackgroundSelection;
class INTupleWriterSvc;
class TTree;
class TFile;
class TLeaf;


/** @class InterleaveAlg
    @brief Gaudi Alg that will fill the ntuple with a background event
 
    It is expected to be invoked for every non-triggering event, on a branch for the Trigger
    algorithm sequence.

 */


class InterleaveAlg : public Algorithm {
public:
    InterleaveAlg(const std::string& name, ISvcLocator* pSvcLocator);
    StatusCode initialize();
    StatusCode execute();
    StatusCode finalize();

    static double triggerRate(){
        return s_triggerRate;
    }

    
private: 
    void copyEventInfo(TTree* t);

    double magneticLatitude();

    BackgroundSelection* m_selector;

    /// access the ntupleWriter service to write out to ROOT ntuples
    INTupleWriterSvc * m_rootTupleSvc;

    StringProperty m_treeName; ///< name of the tree to process
    StringProperty m_rootFile; ///< root file to sample events from
    int m_count;   ///< number of processed events
    int m_downlink;

    TTree * m_meritTuple;

    TLeaf * m_magLatLeaf;


    static double s_triggerRate;


};

#endif
