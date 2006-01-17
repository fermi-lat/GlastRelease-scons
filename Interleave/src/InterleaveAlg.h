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
class ILivetimeSvc;
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

    static double rate(){
        return s_rate;
    }

    
private: 
    void copyEventInfo(); ///< set up 

    double magneticLatitude();

    BackgroundSelection* m_selector;

    /// access the ntupleWriter service to write out to ROOT ntuples
    INTupleWriterSvc * m_rootTupleSvc;

    ILivetimeSvc * m_LivetimeSvc;

    StringProperty m_treeName; ///< name of the tree to process
    StringProperty m_rootFile; ///< root file to sample events from
    StringArrayProperty m_disableList;
    int m_count;   ///< number of processed events
    int m_downlink;

    TTree * m_meritTuple;

    TLeaf * m_magLatLeaf;
    TLeaf * m_runLeaf;
    TLeaf * m_eventLeaf;


    static double s_rate;


};

#endif
