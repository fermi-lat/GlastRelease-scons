/** @file AtwoodTrees.cxx
@brief Implement tree definition and evaluation 

$Header$

*/
#include "GlastClassify/AtwoodTrees.h"
#include "GlastClassify/TreeFactory.h"
#include "TreeAnalysis.h"
#include "src/xmlBuilders/xmlTreeAnalysisFactory.h"

#include "classifier/DecisionTree.h"

#include <sstream>
#include <cassert>

/* 
*/
using namespace GlastClassify;

//_________________________________________________________________________

AtwoodTrees::AtwoodTrees(ITupleInterface& tuple, std::ostream& log, std::string treepath)
                       : m_log(log)
{
    // these are used for preliminary cuts to select the tree to use
    m_TkrNumTracks = tuple.getItem("TkrNumTracks");
    m_CalEnergyRaw = tuple.getItem("CalEnergyRaw");
    m_CalCsIRLn    = tuple.getItem("CalCsIRLn");
    m_EvtEventId   = tuple.getItem("EvtEventId");

    // New items to create or override
    // Add Bill's tuple items so we can start some comparisons
    tuple.addItem("CTBbestEneProb",   m_bestEnergyProb);
    tuple.addItem("CTBprofileProb",   m_profileProb);
    tuple.addItem("CTBlastLayerProb", m_lastLayerProb);
    tuple.addItem("CTBtrackerProb" ,  m_trackerProb);
    tuple.addItem("CTBparamProb",     m_paramProb);
    tuple.addItem("CTBBestEnergy",    m_CTBestEnergy);
    tuple.addItem("CTBdeltaEoE",      m_CTBdeltaEoE);
    tuple.addItem("CTBVTX",           m_VTX);
    tuple.addItem("CTBCORE",          m_CORE);
    tuple.addItem("CTBGAM",           m_GAM);
    
    m_executeTreeCnt = 0;
    m_goodVals       = 0;
    m_caughtVals     = 0;

    //m_xmlFactory = new GlastClassify::xmlTreeFactory(treepath, tuple);
    xmlTreeAnalysisFactory treeFactory(treepath, tuple);

    m_treeAnalysis = treeFactory.buildTreeAnalysis();

    //Testing...
    std::ofstream outFile("IMsheetTest.txt");
    m_treeAnalysis->print(outFile);
    outFile.close();

    return;
}

AtwoodTrees::~AtwoodTrees() 
{
    delete m_treeAnalysis;
}


//_________________________________________________________________________

bool AtwoodTrees::execute()
{
    // initialize CT output variables
    m_bestEnergyProb =  0.;
    m_profileProb    =  0.;
    m_lastLayerProb  =  0.;
    m_trackerProb    =  0.;
    m_paramProb      =  0.;
    m_CTBestEnergy   =  0.;
    m_CTBdeltaEoE    = -2.;
    m_VTX            =  0.;
    m_CORE           =  0.;
    m_GAM            =  0.;

    double tkrNumTracks = *m_TkrNumTracks;
    double calenergy    = *m_CalEnergyRaw;
    double calCsiRln    = *m_CalCsIRLn;
    //double eventId      = *m_EvtEventId;

    // These are the "standard" selection cuts
    if( calenergy <5. || calCsiRln < 4. || tkrNumTracks < 1) return false; 

    m_treeAnalysis->execute();

    m_executeTreeCnt++;

    // Recover the results?
    try
    {
        // Retrieve the energy classification results
        double bestEnergyProb = m_treeAnalysis->getTupleVal("BestEnergyProb");
        double ProfileProb    = m_treeAnalysis->getTupleVal("ProfileProb");
        double lastLayerProb  = m_treeAnalysis->getTupleVal("LastLayerProb");
        double trackerProb    = m_treeAnalysis->getTupleVal("TrackerProb");
        double ParamProb      = m_treeAnalysis->getTupleVal("ParamProb");
        double ctBestEnergy   = m_treeAnalysis->getTupleVal("BestEnergy");
        double ctDeltaEoE     = m_treeAnalysis->getTupleVal("BestDeltaEoE");
        m_bestEnergyProb      = bestEnergyProb;
        m_profileProb         = ProfileProb;
        m_lastLayerProb       = lastLayerProb;
        m_trackerProb         = trackerProb;
        m_paramProb           = ParamProb;
        m_CTBestEnergy        = ctBestEnergy;
        m_CTBdeltaEoE         = ctDeltaEoE;

        // Probability that event was vertexed
        double VTX            = m_treeAnalysis->getTupleVal("VTX");
        m_VTX                 = VTX;

        // Probability that event was "in the core"
        double CORE           = m_treeAnalysis->getTupleVal("CORE");
        m_CORE                = CORE;

        // Background rejection probability
        double GAM            = m_treeAnalysis->getTupleVal("GAM");
        m_GAM                 = GAM;

        m_goodVals++;
    }
    catch(...)
    {
        // Keeps on executing;
        m_caughtVals++;
    }

    double dWriteTupleRow = m_treeAnalysis->getTupleVal("WriteTupleRow");

    bool writeTupleRow = dWriteTupleRow != 0. ? true : false;

    return writeTupleRow;
}
