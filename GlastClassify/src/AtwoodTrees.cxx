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
    tuple.addItem("CTBAcdLowerTileCount",  m_acdLowerTileCount);
    tuple.addItem("CTBAcdUpperTileCount",  m_acdUpperTileCount);
    tuple.addItem("CTBBestPSFerr",         m_bestPsfErr);
    tuple.addItem("CTBBestXDir",           m_bestXDir);
    tuple.addItem("CTBBestYDir",           m_bestYDir);
    tuple.addItem("CTBBestZDir",           m_bestZDir);
    tuple.addItem("CTBBestDeltaEoE",       m_bestDeltaEoE);
    tuple.addItem("CTBBestEnergy",         m_bestEnergy);
    tuple.addItem("CTBBestEnergyProb",     m_bestEnergyProb);
    tuple.addItem("CTBCORE",               m_CORE);
    tuple.addItem("CTBCalFrontBackRatio",  m_calFrontBackRatio);
    tuple.addItem("CTBCalMaxXtalRatio",    m_calMaxXtalRatio);
    tuple.addItem("CTBGAM",                m_GAM);
    tuple.addItem("CTBGoodEnergy",         m_goodEnergy);
    tuple.addItem("CTBLastLayerProb",      m_lastLayerProb);
    tuple.addItem("CTBParamProb",          m_paramProb);
    tuple.addItem("CTBProfileProb",        m_profileProb);
    tuple.addItem("CTBTkrEnergyFrac",      m_tkrEnergyFrac);
    tuple.addItem("CTBTkrLATEdge",         m_tkrLATEdge);
    tuple.addItem("CTBTrackerProb" ,       m_trackerProb);
    tuple.addItem("CTBVTX",                m_VTX);
    
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
    m_acdLowerTileCount   = 0.;
    m_acdUpperTileCount   = 0.;
    m_bestPsfErr          = 0.;
    m_bestXDir            = 0.;
    m_bestYDir            = 0.;
    m_bestZDir            = 0.;
    m_bestDeltaEoE        = 0.;
    m_bestEnergy          = 0.;
    m_bestEnergyProb      = 0.;
    m_CORE                = 0.;
    m_calFrontBackRatio   = 0.;
    m_calMaxXtalRatio     = 0.;
    m_evtLogEnergyRaw     = 0.;
    m_GAM                 = 0.;
    m_goodEnergy          = 0.;
    m_lastLayerProb       = 0.;
    m_paramProb           = 0.;
    m_profileProb         = 0.;
    m_tkrEnergyFrac       = 0.;
    m_tkrLATEdge          = 0.;
    m_trackerProb         = 0.;
    m_VTX                 = 0.;

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
        m_acdLowerTileCount   = m_treeAnalysis->getTupleVal("AcdLowerTileCount");
        m_acdUpperTileCount   = m_treeAnalysis->getTupleVal("AcdUpperTileCount");
        m_bestPsfErr          = m_treeAnalysis->getTupleVal("Best.PSF.err");
        m_bestXDir            = m_treeAnalysis->getTupleVal("Best.XDir");
        m_bestYDir            = m_treeAnalysis->getTupleVal("Best.YDir");
        m_bestZDir            = m_treeAnalysis->getTupleVal("Best.ZDir");
        m_bestDeltaEoE        = m_treeAnalysis->getTupleVal("BestDeltaEoE");
        m_bestEnergy          = m_treeAnalysis->getTupleVal("BestEnergy");
        m_bestEnergyProb      = m_treeAnalysis->getTupleVal("BestEnergyProb");
        m_CORE                = m_treeAnalysis->getTupleVal("CORE");
        m_calFrontBackRatio   = m_treeAnalysis->getTupleVal("CalFrontBackRatio");
        m_calMaxXtalRatio     = m_treeAnalysis->getTupleVal("CalMaxXtalRatio");
        m_evtLogEnergyRaw     = m_treeAnalysis->getTupleVal("EvtLogEnergyRaw");
        m_GAM                 = m_treeAnalysis->getTupleVal("GAM");
        m_goodEnergy          = m_treeAnalysis->getTupleVal("GoodEnergy");
        m_lastLayerProb       = m_treeAnalysis->getTupleVal("LastLayerProb");
        m_paramProb           = m_treeAnalysis->getTupleVal("ParamProb");
        m_profileProb         = m_treeAnalysis->getTupleVal("ProfileProb");
        m_tkrEnergyFrac       = m_treeAnalysis->getTupleVal("TkrEnergyFrac");
        m_tkrLATEdge          = m_treeAnalysis->getTupleVal("TkrLATEdge");
        m_trackerProb         = m_treeAnalysis->getTupleVal("TrackerProb");
        m_VTX                 = m_treeAnalysis->getTupleVal("VTX");

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
