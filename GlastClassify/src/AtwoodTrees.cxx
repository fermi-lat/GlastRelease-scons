/** @file AtwoodTrees.cxx
@brief Implement tree definition and evaluation 

$Header$

*/
#include "GlastClassify/AtwoodTrees.h"
//THB#include "GlastClassify/TreeFactory.h"
#include "GlastClassify/xmlTreeFactory.h"

#include "classifier/DecisionTree.h"

#include <sstream>
#include <cassert>

/* 
*/
using namespace GlastClassify;

namespace {

    // Convenient identifiers used for the nodes
    enum{
        ENERGY_PARAM,     
        ENERGY_LASTLAYER, 
        ENERGY_PROFILE,   
        ENERGY_TRACKER,   
        VERTEX_THIN, VERTEX_THICK,
        PSF_VERTEX_THIN, 
        PSF_VERTEX_THICK,
        PSF_TRACK_THIN,  
        PSF_TRACK_THICK, 
        GAMMA_VERTEX_HIGH, GAMMA_VERTEX_MED, GAMMA_VERTEX_THIN, GAMMA_VERTEX_THICK,
        GAMMA_TRACK_HIGH,  GAMMA_TRACK_MED,  GAMMA_TRACK_THIN,  GAMMA_TRACK_THICK,
        NODE_COUNT, // stop here
    };

    /** table to correlate indeces with 
      */
    //___________________________________________________________________________

    class CTinfo { 
    public:
        int id;           // unique ID for local identification
        std::string name; // the name of the DecisionTree
    };
    // these have to correspond to the folder names
    CTinfo imNodeInfo[] = {
        { ENERGY_PARAM,      "energy/param" },
        { ENERGY_LASTLAYER,  "energy/lastlayer" },
        { ENERGY_PROFILE,    "energy/profile" },
        { ENERGY_TRACKER,    "energy/tracker" },

        { VERTEX_THIN,       "vertex/thin"},
        { VERTEX_THICK,      "vertex/thick"},

        { PSF_VERTEX_THIN,   "psf/vertex/thin"}, 
        { PSF_VERTEX_THICK,  "psf/vertex/thick"},
        { PSF_TRACK_THIN,    "psf/track/thin"},
        { PSF_TRACK_THICK,   "psf/track/thick"},
        { GAMMA_VERTEX_HIGH, "gamma/vertex/highcal"},
        { GAMMA_VERTEX_MED,  "gamma/vertex/medcal"},
        { GAMMA_VERTEX_THIN, "gamma/vertex/thin"},
        { GAMMA_VERTEX_THICK,"gamma/vertex/thick"},
        { GAMMA_TRACK_HIGH,  "gamma/track/highcal"},
        { GAMMA_TRACK_MED,   "gamma/track/medcal"},
        { GAMMA_TRACK_THIN,  "gamma/track/thin"},
        { GAMMA_TRACK_THICK, "gamma/track/thick"},
    };


}  // anonymous namespace

//_________________________________________________________________________

AtwoodTrees::AtwoodTrees( 
     ITupleInterface& tuple, 
     std::ostream&             log,       
     std::string               treepath
     )
     : m_log(log)
{
    // these are used for preliminary cuts to select the tree to use
    m_Tkr1FirstLayer = tuple.getItem("Tkr1FirstLayer");
    m_CalEnergyRaw=tuple.getItem("CalEnergyRaw"     );
    m_CalTotRLn   =tuple.getItem("CalTotRLn"        );
    m_VtxAngle    =tuple.getItem("VtxAngle"         ); 

    // the energy estimates
    m_EvtEnergyCorr = tuple.getItem("EvtEnergyCorr"  );
    m_CalCfpEnergy=   tuple.getItem("CalCfpEnergy"   );
    m_CalLllEnergy =  tuple.getItem("CalLllEnergy"   );
    m_CalTklEnergy=   tuple.getItem("CalTklEnergy"   );



    // New items to create or override
    // create new float TupleItem objects: will be automatically added to the overall tuple
    tuple.addItem("CTgoodCal",  m_goodCalProb);
    tuple.addItem("CTvertex",   m_vtxProb);
    tuple.addItem("CTgoodPsf",  m_goodPsfProb);
    tuple.addItem("CTgamma" ,   m_gammaProb);
    tuple.addItem("CTgammaType",m_gammaType);
    tuple.addItem("BestEnergy", m_BestEnergy);

    //m_factory = new GlastClassify::TreeFactory(treepath, tuple);
    m_factory = new GlastClassify::xmlTreeFactory(treepath, tuple);

    for( unsigned int i=0; i<NODE_COUNT; ++i)
    {
        const ITreeFactory::ITree& tree = (*m_factory)(imNodeInfo[i].name);

        //std::string sOutFileRoot = treepath + "/" + imNodeInfo[i].name;

        //dynamic_cast<const xmlTreeFactory::Tree&>(tree).printFile(sOutFileRoot);
    }
}

    //_________________________________________________________________________ 

    bool AtwoodTrees::useVertex()const
    {
        return *m_VtxAngle>0 && m_vtxProb >0.5;
    }

    //_________________________________________________________________________

    void AtwoodTrees::execute()
    {

        // initialize CT output variables
        m_goodPsfProb=0;
        m_vtxProb  = m_gammaProb = m_goodCalProb= 0;
        m_gammaType = -1;

        double calenergy = *m_CalEnergyRaw;
        if( calenergy <5. || *m_CalTotRLn < 4.) return; // the "NoCAL" case that we cannot deal with

        double energymeasure[] = {*m_EvtEnergyCorr, *m_CalCfpEnergy, *m_CalLllEnergy, *m_CalTklEnergy};
        int ctree_index[] = {ENERGY_PARAM, ENERGY_PROFILE, ENERGY_LASTLAYER, ENERGY_TRACKER};
        double bestprob(0);
        // ============> wire in only param here (for now) <===================
        for( int i =0; i<1; ++i){
            double prob = energymeasure[i]>0? m_factory->evaluate(ctree_index[i]) : 0;
            if( prob>bestprob){
                bestprob = prob;
                m_BestEnergy = energymeasure[i];
            }
        }

        // assign tuple items
        m_goodCalProb = bestprob;

        // evaluate the rest only if cal prob ok (should this be wired in??? (Bill now has 0.1)
        if( m_goodCalProb<0.25 )   return;

        // selection of energy range for gamma trees
        enum {CAL_LOW, CAL_MED, CAL_HIGH};
        int cal_type;
        if(     calenergy <  350) cal_type = CAL_LOW;
        else if(calenergy < 3500) cal_type = CAL_MED;
        else                      cal_type = CAL_HIGH;


        // select vertex-vs-track, and corresponding trees for good psf, background rejection

        int psfType=0, gammaType=0;
        if( *m_Tkr1FirstLayer > 5 ) { // thin

            m_vtxProb = m_factory->evaluate(VERTEX_THIN); 
            if( useVertex() ) {
                psfType = PSF_VERTEX_THIN ;
                if(      cal_type==CAL_HIGH) gammaType=GAMMA_VERTEX_HIGH;
                else if( cal_type==CAL_MED)  gammaType=GAMMA_VERTEX_MED;
                else                         gammaType=GAMMA_VERTEX_THIN;

            }else{ //track
                psfType = PSF_TRACK_THIN;
                if(      cal_type==CAL_HIGH) gammaType=GAMMA_TRACK_HIGH;
                else if( cal_type==CAL_MED)  gammaType=GAMMA_TRACK_MED;
                else                         gammaType=GAMMA_TRACK_THIN;
            }
        }else { //thick
            m_vtxProb = m_factory->evaluate(VERTEX_THICK);
           if( useVertex() ) {
                psfType = PSF_VERTEX_THICK ;
                if(      cal_type==CAL_HIGH) gammaType=GAMMA_VERTEX_HIGH;
                else if( cal_type==CAL_MED)  gammaType=GAMMA_VERTEX_MED;
                else                         gammaType=GAMMA_VERTEX_THICK;

            }else{ //track
                psfType = PSF_TRACK_THICK;
                if(      cal_type==CAL_HIGH) gammaType=GAMMA_TRACK_HIGH;
                else if( cal_type==CAL_MED)  gammaType=GAMMA_TRACK_MED;
                else                         gammaType=GAMMA_TRACK_THICK;
            }
        }

        // now evalute the appropriate trees
        m_goodPsfProb = m_factory->evaluate(psfType);

        m_gammaProb   = m_factory->evaluate(gammaType);
        m_gammaType   = gammaType-GAMMA_VERTEX_HIGH; // will be 0-7
    }

    //_________________________________________________________________________

    AtwoodTrees::~AtwoodTrees()
    {
        delete m_factory;
    }
