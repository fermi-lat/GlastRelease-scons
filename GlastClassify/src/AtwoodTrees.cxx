/** @file AtwoodTrees.cxx
@brief Implement tree definition and evaluation 

$Header$

*/
#include "AtwoodTrees.h"
#include "TreeAnalysis.h"
#include "xmlBuilders/xmlTreeAnalysisFactory.h"

#include "classifier/DecisionTree.h"
#include "facilities/Util.h" // for expandEnvVar    


#include <sstream>
#include <cassert>

namespace {
    // This should be changed each time there is a new file, to make it defautl
    std::string default_xml("$(GLASTCLASSIFYROOT)/xml/Pass5_Analysis_Complete_4_PSep.xml");
}
/* 
*/
using namespace GlastClassify;

//_________________________________________________________________________

AtwoodTrees::AtwoodTrees(ITupleInterface& tuple, std::ostream& log, std::string imfile)
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

    m_eventId          = tuple.getItem("EvtEventId");
    m_run              = tuple.getItem("EvtRun");

    // Default values as of 4/5/2007 at request of Bill Atwood
    m_calEnergyCut     = 5.;
    m_csiRadLenCut     = 4.;
    m_numTracksCut     = 0.;
    
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

    // Unfortunately, the following is no longer the way "CTB" variables get added to the 
    // output ntuple. At the request of Bill Atwood, ALL variables in the IM analysis which
    // are prefixed by "CTB" are now output to the ntuple. Unfortunately, I have no idea 
    // how we are going to document this... The list below is more or less what gets written
    // out so perhaps it will still be useful...
    // New items to create or override
    // Add Bill's tuple items so we can start some comparisons
    /*
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
    */
    
    m_executeTreeCnt = 0;
    m_goodVals       = 0;
    m_caughtVals     = 0;

    //m_xmlFactory = new GlastClassify::xmlTreeFactory(treepath, tuple);
    if( imfile.empty() ){
        imfile = default_xml; 
    }

    facilities::Util::expandEnvVar(&imfile);

    log << "GlastClassify::AtwoodTrees-- Loading Atwood's IM classification trees from " << imfile ;

    xmlTreeAnalysisFactory treeFactory(imfile, tuple);

    m_treeAnalysis = treeFactory.buildTreeAnalysis();

    log << "\n\t\tloaded " << treeFactory.nodeCount() <<" nodes";

    //Testing...
    //std::ofstream outFile("IMsheetTest.txt");
    //m_treeAnalysis->print(outFile);
    //outFile.close();

    return;
}

AtwoodTrees::~AtwoodTrees() 
{
    delete m_treeAnalysis;
}


//_________________________________________________________________________

bool AtwoodTrees::execute()
{
    // Output flag default is to not write row (if pruning)
    bool writeTupleRow = false;

    // initialize CT output variables
    m_bestEnergyProb    = 0.;
    m_CORE              = 0.;
    m_evtLogEnergyRaw   = 0.;

    double tkrNumTracks = *m_TkrNumTracks;
    double calenergy    = *m_CalEnergyRaw;
    double calCsiRln    = *m_CalCsIRLn;

    int    eventId      = *m_eventId;
    int    run          = *m_run;

    // Always zero the CTB output values in case cuts below fail
    m_treeAnalysis->zeroCTvals();

    // These are the "standard" selection cuts
    if( calenergy > m_calEnergyCut && calCsiRln > m_csiRadLenCut && tkrNumTracks > m_numTracksCut)
    {
        m_treeAnalysis->execute();

        m_executeTreeCnt++;

        // Recover the results?
        try
        {
            // Retrieve the energy classification results (needed below)
            m_bestEnergyProb  = m_treeAnalysis->getTupleVal("CTBBestEnergyProb");
            m_CORE            = m_treeAnalysis->getTupleVal("CTBCORE");
            m_evtLogEnergyRaw = m_treeAnalysis->getTupleVal("EvtLogEnergyRaw");

            m_goodVals++;
        }
        catch(...)
        {
            // Keeps on executing;
            m_caughtVals++;
        }

        double FilterStatus_HI = *m_FilterStatus_HI;

        FilterStatus_HI = 0;

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

                if (!(AcdCornerDoca > -5 && AcdCornerDoca < 50))
                {
                    // Finally, check the result of running the Analysis Sheet
                    double dWriteTupleRow = m_treeAnalysis->getTupleVal("WriteTupleRow");

                    if (dWriteTupleRow != 0.) writeTupleRow = true;
                }
            }
        }
    }

    // Last step, always copy CTB back to ntuple
    m_treeAnalysis->storeCTvals();

    return writeTupleRow;
}
