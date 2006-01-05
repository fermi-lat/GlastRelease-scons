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
    m_TkrNumTracks    = tuple.getItem("TkrNumTracks");
    m_CalEnergyRaw    = tuple.getItem("CalEnergyRaw");
    m_CalCsIRLn       = tuple.getItem("CalCsIRLn");
    m_EvtEventId      = tuple.getItem("EvtEventId");
    m_FilterStatus_HI = tuple.getItem("FilterStatus_HI");

    m_AcdActiveDist3D  = tuple.getItem("AcdActiveDist3D");
    m_AcdRibbonActDist = tuple.getItem("AcdRibbonActDist");
    m_AcdCornerDoca    = tuple.getItem("AcdCornerDoca");
    m_Tkr1SSDVeto      = tuple.getItem("Tkr1SSDVeto");

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
    tuple.addItem("CTBCalDocaAngle",       m_calDocaAngle);
    tuple.addItem("CTBCalFrontBackRatio",  m_calFrontBackRatio);
    tuple.addItem("CTBCalMaxXtalRatio",    m_calMaxXtalRatio);
    tuple.addItem("CTBCalTransTCCD",       m_calTransTCCD);
    tuple.addItem("CTBGAM",                m_GAM);
    //tuple.addItem("CTBGoodEnergy",         m_goodEnergy); **REMOVE**
    tuple.addItem("CTBLastLayerProb",      m_lastLayerProb);
    tuple.addItem("CTBParamProb",          m_paramProb);
    tuple.addItem("CTBProfileProb",        m_profileProb);
    tuple.addItem("CTBTkrCoreCalDoca",     m_tkrCoreCalDoca);
    tuple.addItem("CTBTkrEnergyFrac",      m_tkrEnergyFrac);
    tuple.addItem("CTBTkrLATEdge",         m_tkrLATEdge);
    tuple.addItem("CTBTkrSHRCalAngle" ,    m_tkrSHRCalAngle);
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
    m_calDocaAngle        = 0.;        // Added 1/3/06
    m_calFrontBackRatio   = 0.;
    m_calMaxXtalRatio     = 0.;
    m_calTransTCCD        = 0.;        // Added 1/3/06
    m_evtLogEnergyRaw     = 0.;
    m_GAM                 = 0.;
    //m_goodEnergy          = 0.;
    m_lastLayerProb       = 0.;
    m_paramProb           = 0.;
    m_profileProb         = 0.;
    m_tkrCoreCalDoca      = 0.;      // Added 1/3/06
    m_tkrEnergyFrac       = 0.;
    m_tkrSHRCalAngle      = 0.;      // Added 1/3/06
    m_tkrLATEdge          = 0.;
    m_trackerProb         = 0.;
    m_VTX                 = 0.;

    double tkrNumTracks = *m_TkrNumTracks;
    double calenergy    = *m_CalEnergyRaw;
    double calCsiRln    = *m_CalCsIRLn;
    //double eventId      = *m_EvtEventId;

    // These are the "standard" selection cuts
    if( calenergy <= 5. || calCsiRln <= 4. || tkrNumTracks < 1) return false; 

    m_treeAnalysis->execute();

    m_executeTreeCnt++;

    // Recover the results?
    try
    {
        // Retrieve the energy classification results
        m_acdLowerTileCount   = m_treeAnalysis->getTupleVal("CTBAcdLowerTileCount");
        m_acdUpperTileCount   = m_treeAnalysis->getTupleVal("CTBAcdUpperTileCount");
        m_bestPsfErr          = m_treeAnalysis->getTupleVal("Best.PSF.err");
        m_bestXDir            = m_treeAnalysis->getTupleVal("Best.XDir");
        m_bestYDir            = m_treeAnalysis->getTupleVal("Best.YDir");
        m_bestZDir            = m_treeAnalysis->getTupleVal("Best.ZDir");
        m_bestDeltaEoE        = m_treeAnalysis->getTupleVal("BestDeltaEoE");
        m_bestEnergy          = m_treeAnalysis->getTupleVal("BestEnergy");
        m_bestEnergyProb      = m_treeAnalysis->getTupleVal("BestEnergyProb");
        m_CORE                = m_treeAnalysis->getTupleVal("CORE");
        m_calDocaAngle        = m_treeAnalysis->getTupleVal("CTBCalDocaAngle");        // Added 1/3/06
        m_calFrontBackRatio   = m_treeAnalysis->getTupleVal("CTBCalFrontBackRatio");
        m_calMaxXtalRatio     = m_treeAnalysis->getTupleVal("CTBCalMaxXtalRatio");
        m_calTransTCCD        = m_treeAnalysis->getTupleVal("CTBCalTransTCCD");        // Added 1/3/06
        m_evtLogEnergyRaw     = m_treeAnalysis->getTupleVal("EvtLogEnergyRaw");
        m_GAM                 = m_treeAnalysis->getTupleVal("GAM");
//        m_goodEnergy          = m_treeAnalysis->getTupleVal("GoodEnergy");
        m_lastLayerProb       = m_treeAnalysis->getTupleVal("LastLayerProb");
        m_paramProb           = m_treeAnalysis->getTupleVal("ParamProb");
        m_profileProb         = m_treeAnalysis->getTupleVal("ProfileProb");
        m_tkrCoreCalDoca      = m_treeAnalysis->getTupleVal("CTBTkrCoreCalDoca");      // Added 1/3/06
        m_tkrEnergyFrac       = m_treeAnalysis->getTupleVal("CTBTkrEnergyFrac");
        m_tkrLATEdge          = m_treeAnalysis->getTupleVal("CTBTkrLATEdge");
        m_tkrSHRCalAngle      = m_treeAnalysis->getTupleVal("CTBTkrSHRCalAngle");      // Added 1/3/06
        m_trackerProb         = m_treeAnalysis->getTupleVal("TrackerProb");
        m_VTX                 = m_treeAnalysis->getTupleVal("VTX");

        m_goodVals++;
    }
    catch(...)
    {
        // Keeps on executing;
        m_caughtVals++;
    }


    // Cuts for Special Bill Run
    bool writeTupleRow = false;

    double FilterStatus_HI = *m_FilterStatus_HI;

    // First cuts on Filter status and failures for energy and tails
    if (FilterStatus_HI == 0 && m_bestEnergyProb > 0.1 && m_CORE > 0.1)
    {
        double AcdActiveDist3D  = *m_AcdActiveDist3D;
        double AcdRibbonActDist = *m_AcdRibbonActDist;
        double Tkr1SSDVeto      = *m_Tkr1SSDVeto;

        // A series of selections on the ACD 
        if (!((AcdActiveDist3D > 0 || AcdRibbonActDist > 0) && Tkr1SSDVeto < 2))
        {
            double AcdCornerDoca = *m_AcdCornerDoca;

            if (!(AcdCornerDoca > -5 && AcdCornerDoca < 50 && m_tkrLATEdge < 100))
            {
                // Finally, check the result of running the Analysis Sheet
                double dWriteTupleRow = m_treeAnalysis->getTupleVal("WriteTupleRow");

                if (dWriteTupleRow != 0.) writeTupleRow = true;
            }
        }
    }

    return writeTupleRow;
}
