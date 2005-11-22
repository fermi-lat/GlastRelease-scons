/** @file AtwoodLikeTrees.cxx
@brief Implement tree definition and evaluation 

$Header$

*/
#include "GlastClassify/AtwoodLikeTrees.h"

#include "classifier/DecisionTree.h"

#include <sstream>
#include <cassert>
#include <map>
#include <cmath>

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

    template< typename T>
    class LocalValue : public GlastClassify::Item {
    public:
        LocalValue(const T & val):m_val(val){}
        operator double()const{return m_val;}
    private:
        const T& m_val;
    };

    inline double sqr(double x){return x*x;}

}  // anonymous namespace

//_________________________________________________________________________

AtwoodLikeTrees::AtwoodLikeTrees( 
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
    
        // for evaluation of local variables
    m_TkrEnergyCorr = tuple.getItem("TkrEnergyCorr"  );
    m_Tkr1Hits      = tuple.getItem("Tkr1Hits");
    m_Tkr2Hits      = tuple.getItem("Tkr1Hits");
    m_Tkr1X0        = tuple.getItem("Tkr1X0");
    m_Tkr1Y0        = tuple.getItem("Tkr1Y0");
    m_TkrTotalHits  = tuple.getItem("TkrTotalHits");
    m_CalXtalMaxEne = tuple.getItem("CalXtalMaxEne");
    

    // New items to create or override
    // create new float TupleItem objects: will be automatically added to the overall tuple
#if 0
    std::string prefix("IM"); 
#else
    std::string prefix("CT"); 
#endif
    tuple.addItem(prefix+"goodCal",  m_goodCalProb);
    tuple.addItem(prefix+"vertex",   m_vtxProb);
    tuple.addItem(prefix+"goodPsf",  m_goodPsfProb);
    tuple.addItem(prefix+"gamma" ,   m_gammaProb);
    tuple.addItem(prefix+"gammaType",m_gammaType);
    tuple.addItem(prefix+"BestEnergy", m_BestEnergy);

    // set up list of locally defined  variables for access by the trees
                    
    std::map<std::string, const Item*> localvals;
    localvals["EvtLogEnergyRaw"] = new LocalValue<float>(m_EvtLogEnergyRaw);
    localvals["BestLogEnergy"]   = new LocalValue<float>(m_BestLogEnergy);
    localvals["TkrEnergyFrac"]   = new LocalValue<float>(m_TkrEnergyFrac);
    localvals["TkrTotalHitsNorm"]= new LocalValue<float>(m_TkrTotalHitsNorm);
    localvals["TkrTotalHitRatio"]= new LocalValue<float>(m_TkrTotalHitRatio);
    localvals["BestEnergyProb"]  = new LocalValue<float>(m_BestEnergyProb);
    localvals["Tkr12DiffHits"]   = new LocalValue<float>(m_Tkr12DiffHits);
    localvals["PSFEneProbPrd"]   = new LocalValue<float>(m_PSFEneProbPrd);
    localvals["CalMaxXtalRatio"] = new LocalValue<float>(m_CalMaxXtalRatio);
    localvals["TkrLATEdge"]      = new LocalValue<float>(m_TkrLATEdge);

    // these are aliases -- Tracy has the following maping
    /*
    m_ImToTobyOutMap["Pr(GoodEnergy)"] = "CTgoodCal";
    m_ImToTobyOutMap["Pr(CORE)"]       = "CTgoodPsf";
    m_ImToTobyOutMap["Pr(VTX)"]        = "CTvertex";
    m_ImToTobyOutMap["Pr(GAM)"]        = "CTgamma";
    */

    localvals["CTvertex"]        = new LocalValue<double>(m_vtxProb);
    localvals["CTgoodPsf"]       = new LocalValue<double>(m_goodPsfProb);
    // only seem to need these

    m_factory = new GlastClassify::TreeFactory(treepath, tuple, localvals);

    for( unsigned int i=0; i<NODE_COUNT; ++i){
        (*m_factory)(imNodeInfo[i].name);
    }
}

    //_________________________________________________________________________ 

    bool AtwoodLikeTrees::useVertex()const
    {
        return *m_VtxAngle>0 && m_vtxProb >0.5;
    }

    //_________________________________________________________________________

    void AtwoodLikeTrees::execute()
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
        for( int i =0; i<4; ++i){
            double prob = energymeasure[i]>0? m_factory->evaluate(ctree_index[i]) : 0;
            if( prob>bestprob){
                bestprob = prob;
                m_BestEnergy = energymeasure[i];
            }
        }
//        BestLogEnergy = log10(m_BestEnergy);

        // assign tuple items
        m_goodCalProb = bestprob;

        // evaluate the rest only if cal prob ok (should this be wired in??? (Bill now has 0.1)
        if( m_goodCalProb<0.25 )   return;

        m_BestLogEnergy = log10(m_BestEnergy);
        m_BestEnergyProb = bestprob;
        m_TkrEnergyFrac = *m_TkrEnergyCorr/ *m_EvtEnergyCorr;
        m_TkrTotalHitsNorm = *m_TkrTotalHits/(*m_Tkr1FirstLayer-1);
        m_TkrTotalHitRatio = (*m_Tkr1Hits+*m_Tkr2Hits) / (*m_TkrTotalHits);
        m_Tkr12DiffHits = *m_Tkr1Hits - *m_Tkr2Hits;
        m_TkrLATEdge = 740. - std::max(fabs(*m_Tkr1X0), fabs(*m_Tkr1Y0));
        m_CalMaxXtalRatio = *m_CalXtalMaxEne / *m_CalEnergyRaw;



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

        m_PSFEneProbPrd = sqrt(sqr(m_BestEnergyProb)+sqr(m_goodPsfProb));

        m_gammaProb   = m_factory->evaluate(gammaType);
        m_gammaType   = gammaType-GAMMA_VERTEX_HIGH; // will be 0-7
    }

    //_________________________________________________________________________

    AtwoodLikeTrees::~AtwoodLikeTrees()
    {
        delete m_factory;
    }
