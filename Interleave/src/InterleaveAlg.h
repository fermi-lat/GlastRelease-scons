/** @file InterleaveAlgcxx

@brief declaration  of the class InterleaveAlg

$Header$

*/
#ifndef InterleaveAlg_h
#define InterleaveAlg_h


#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/Property.h"
class INTupleWriterSvc;
class TTree;
class TFile;

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
    void copyTreeData(TTree* in, TTree* out);

    /// access the ntupleWriter service to write out to ROOT ntuples
    INTupleWriterSvc * m_rootTupleSvc;

    StringProperty m_treeName; ///< name of the tree to process
    StringProperty m_rootFile; ///< root file to sample events from
    int m_count;   ///< number of processed events
    int m_downlink; ///< number of down linked events
    int m_event; ///< current event number in root file

    TFile * m_tfile;
    TTree * m_backgroundTuple;
    TTree * m_meritTuple;


    static double s_triggerRate;
    static double s_downlinkRate;


};

#endif
