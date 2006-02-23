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
    m_FilterStatus_HI = tuple.getItem("FilterStatus_HI");

    m_AcdActiveDist3D  = tuple.getItem("AcdActiveDist3D");
    m_AcdRibbonActDist = tuple.getItem("AcdRibbonActDist");
    m_AcdCornerDoca    = tuple.getItem("AcdCornerDoca");
    m_Tkr1SSDVeto      = tuple.getItem("Tkr1SSDVeto");
    
/** @page merittuple
@section ctbVars CTB Variables
<table>
    <tr><th> Variable <th> Type  <th> Description				                 
    <tr><td> CTBAcdLowerTileCount
    <td>F<td>  AcdNoSideRow3
    <tr><td> CTBAcdUpperTileCount
    <td>F<td>  AcdNoTop+AcdNoSideRow0+AcdNoSideRow1+AcdSideRow2    
    <tr><td> CTBBestPSFerr
    <td>F<td>  The angle between the best direction and the MC direction (radians)
    <tr><td> CTBBest[X/Y/Z]Dir
    <td>F<td>  Best direction, chosen from vertex and track-1 Solutions (3 direction cosines)  
    <tr><td> CTBBestDeltaEoE
    <td>F<td>  Best Energy Error relative to MC energy D(E)/E    
    <tr><td> CTBBestEnergy
    <td>F<td>  Best Estimated energy from among the 4 methods    
    <tr><td> CTBBestEnergyProb
    <td>F<td>  Energy probability knob.
               Tunes the energy resolution for the selected method    
    <tr><td> CTBBestLogEnergy
    <td>F<td>  Log(CTBBestEnergy) – base 10    
    <tr><td> CTBCORE
    <td>F<td>  Image probability knob. Tunes the image resolution  
    <tr><td> CTBCalDocaAngle
    <td>F<td>  CalTrackDoca + 80*CalTrackAngle    
    <tr><td> CTBCalFrontBackRatio
    <td>F<td>  Probably the ratio of the energies in the front and back halves of the CAL.
    <tr><td> CTBCalMaxXtalRatio
    <td>F<td>  CalXtalMaxEne/CalEnergyRaw    
    <tr><td> CTBCalTransTCCD
    <td>F<td>  CalTransRms + 0.1*(CalTrackDoca - 2.5*Tkr1CoreHC)    
    <tr><td> CTBGAM
    <td>F<td>  Background rejection probability knob. Tunes the background contamination
    <tr><td> CTBLastLayerProb
    <td>F<td>  Probibility for the "corrections" of LastLayer energy method 
               against a fixed functional standard   
    <tr><td> CTBParamProb
    <td>F<td>  Same for Parameter energy method  
    <tr><td> CTBProfileProb
    <td>F<td>  Same for Profile energy method 
    <tr><td> CTBTrackerProb
    <td>F<td>  Same for Tracker energy method
    <tr><td> CTBTkrCoreCalDoca
    <td>F<td>  CalTrackDoca - 2.5*Tkr1CoreHC. Ad hoc background rejection variable    
    <tr><td> CTBTkrEnergyFrac
    <td>F<td>  TkrEnergyCorr/EvtEnergyCorr –  Ad hoc background rejection variable       
    <tr><td> CTBTkrLATEdge
    <td>F<td>  742. - max(abs(Tkr1X0) , abs(Tkr1Y0)). 
               Fiducial Volume Variable, more accurately calculated as Tkr1LATEdge
    <tr><td> CTBTkrSHRCalAngle
    <td>F<td>  CalTrackAngle - 0.2*TkrSurplusHitRatio. Ad hoc background rejection variable
    <tr><td> CTBVTX
    <td>F<td>  Internal probability to select between track-1 and vertex
</table>
*/



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
    tuple.addItem("CTBBestLogEnergy",      m_bestLogEnergy);
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
    m_bestLogEnergy       = 0.;
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
        m_bestDeltaEoE        = m_treeAnalysis->getTupleVal("CTBBestDeltaEoE");
        m_bestEnergy          = m_treeAnalysis->getTupleVal("CTBBestEnergy");
        m_bestEnergyProb      = m_treeAnalysis->getTupleVal("CTBBestEnergyProb");
        m_bestLogEnergy       = m_treeAnalysis->getTupleVal("CTBBestLogEnergy");
        m_bestPsfErr          = m_treeAnalysis->getTupleVal("CTBBestPSFerr");
        m_bestXDir            = m_treeAnalysis->getTupleVal("CTBBestXDir");
        m_bestYDir            = m_treeAnalysis->getTupleVal("CTBBestYDir");
        m_bestZDir            = m_treeAnalysis->getTupleVal("CTBBestZDir");
        m_CORE                = m_treeAnalysis->getTupleVal("CTBCORE");
        m_calDocaAngle        = m_treeAnalysis->getTupleVal("CTBCalDocaAngle");        // Added 1/3/06
        m_calFrontBackRatio   = m_treeAnalysis->getTupleVal("CTBCalFrontBackRatio");
        m_calMaxXtalRatio     = m_treeAnalysis->getTupleVal("CTBCalMaxXtalRatio");
        m_calTransTCCD        = m_treeAnalysis->getTupleVal("CTBCalTransTCCD");        // Added 1/3/06
        m_evtLogEnergyRaw     = m_treeAnalysis->getTupleVal("EvtLogEnergyRaw");
        m_GAM                 = m_treeAnalysis->getTupleVal("CTBGAM");
//        m_goodEnergy          = m_treeAnalysis->getTupleVal("GoodEnergy");
        m_lastLayerProb       = m_treeAnalysis->getTupleVal("LastLayerProb");
        m_paramProb           = m_treeAnalysis->getTupleVal("ParamProb");
        m_profileProb         = m_treeAnalysis->getTupleVal("ProfileProb");
        m_tkrCoreCalDoca      = m_treeAnalysis->getTupleVal("CTBTkrCoreCalDoca");      // Added 1/3/06
        m_tkrEnergyFrac       = m_treeAnalysis->getTupleVal("CTBTkrEnergyFrac");
        m_tkrLATEdge          = m_treeAnalysis->getTupleVal("CTBTkrLATEdge");
        m_tkrSHRCalAngle      = m_treeAnalysis->getTupleVal("CTBTkrSHRCalAngle");      // Added 1/3/06
        m_trackerProb         = m_treeAnalysis->getTupleVal("TrackerProb");
        m_VTX                 = m_treeAnalysis->getTupleVal("CTBVTX");

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
