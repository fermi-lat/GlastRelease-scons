/**@file xmlTreeFactory.cxx

@brief implementation of class xmlTreeFactory

$Header$
*/

#include "GlastClassify/xmlTreeFactory.h"
#include "classifier/DecisionTree.h"
#include "ImActivityNodes/PredictEngineNode.h"
#include "xmlBuilders/ImSheetBuilder.h"
#include <fstream>
#include <cassert>
#include <stdexcept>
#include <algorithm>
#include <cmath>  // for M_PI, among others

#include "xmlBase/XmlParser.h"
#include <xercesc/dom/DOMElement.hpp>
#include "xmlBase/Dom.h"

using namespace GlastClassify;
XERCES_CPP_NAMESPACE_USE


/** @class GleamValues
@brief local definition of class which handles pointers to values
       NOTE: this needs modification - separate out into own class
       (in progress for next check in)
*/
class xmlTreeFactory::GleamValues : public DecisionTree::Values 
{
public:
    GleamValues(const ImSheetBuilder::StringList& names, 
                const xmlTreeFactory::LocalTupleValues& localVals,
                const std::map<std::string,std::string>& imToTobyMap,
                ITreeFactory::ILookupData& lookup) :
                m_localVals(localVals)
    {
        // Loop over list of variable names given to us in names
        //for( DecisionTreeBuilder::StringList::const_iterator it = names.begin();
        for( ImSheetBuilder::StringList::const_iterator it = names.begin();
            it != names.end(); ++it)
        {
            std::string varName = *it;

            // Does this variable need to have its name changed?
            if (imToTobyMap.find(varName) != imToTobyMap.end())
            {
                varName = imToTobyMap.find(varName)->second;
            }

            // Use a "try" to catch the exception if not in the tuple
            try
            {
                std::pair<bool, const void*> entry = lookup(varName);
                if( entry.second == 0 ){
                    throw std::invalid_argument("xmlTreeFactory: did not find variable "+varName);
                }

                m_varAction.push_back(std::pair<std::string, bool>(varName, true));
                m_pval.push_back( entry);
            }
            catch (std::runtime_error&)
            {
                // Is this value a valid variable to be computed?
                if (localVals.isValue(varName))
                {
                    //Need a place holder for indexing to work
                    std::pair<bool, const void*> entry(false, (const void*)0);
                    
                    m_varAction.push_back(std::pair<std::string, bool>(varName, false));
                    m_pval.push_back(entry);
                }
                // Otherwise, an error
                else
                {
                    int j = 0;
                    //throw std::invalid_argument("xmlTreeFactory: did not find variable "+*it);
                }
            }
        }

    }
    /// @brief callback from tree evaluation

    double operator[](int index)const
    {
        double val = 0;
        const std::pair<std::string, bool>& varActionPair = m_varAction[index];

        if (varActionPair.second)
        {
            std::pair<bool, const void*> entry = m_pval[index];
            // now dereference either as a float or a double
            val = entry.first? *(const float*)entry.second : *(const double*)entry.second;
        }
        else
        {
            val = m_localVals.getValue(varActionPair.first);
        }

        return val;
    }
//private:
    std::vector<std::pair<std::string, bool> > m_varAction;
    std::vector<std::pair<bool, const void*> > m_pval;
    const xmlTreeFactory::LocalTupleValues&    m_localVals;
};


xmlTreeFactory::xmlTreeFactory(const std::string& path, ITreeFactory::ILookupData& lookup)
                               : m_lookup(lookup), m_localVals(lookup)
{
    std::string sFileName = path+"/"+"DC2_Analysis.imw";
    
    xmlBase::XmlParser xmlParser;

    m_domDocument = xmlParser.parse(sFileName.c_str());
    
    if(m_domDocument == 0)
    {
        // Error checking usually for a missing file.
        // Or sometimes due to a bad XML document.
        // Remember <P><BR></P> is only valid
        // in poorly formed html, not xml.
        // When we get a schema for UserLibrary, we'll
        // be able to use that as a validator also.
        std::invalid_argument("xmlTreeFactory: could no open input file: " + sFileName);
    }

    // Testing
    m_imSheet = new ImSheetBuilder(m_domDocument);

    //Testing...
    std::ofstream outFile("IMsheetTest.txt");
    m_imSheet->print(outFile);
    outFile.close();

    // Make a map between the names Toby has for the trees, and Bill's names...
    m_TobyToBillMap.clear();

    m_TobyToBillMap["energy/param"]         = "CT: Energy Param";
    m_TobyToBillMap["energy/lastlayer"]     = "CT: Energy Last Layer";
    m_TobyToBillMap["energy/profile"]       = "CT: Energy Profile";
    m_TobyToBillMap["energy/tracker"]       = "CT: Energy Tracker";

    m_TobyToBillMap["vertex/thin"]          = "CT: PSF VTX/1Tkr Thin";
    m_TobyToBillMap["vertex/thick"]         = "CT: PSF VTX/1Tkr Thick";

    m_TobyToBillMap["psf/vertex/thin"]      = "CT: PSF CORE VTX Thin";  
    m_TobyToBillMap["psf/vertex/thick"]     = "CT: PSF CORE VTX Thick"; 
    m_TobyToBillMap["psf/track/thin"]       = "CT: PSF CORE 1TKR Thin"; 
    m_TobyToBillMap["psf/track/thick"]      = "CT:  PSF CORE 1TKR Thick";

    // Don't know what to do here but need at least a placeholder...
    m_TobyToBillMap["gamma/vertex/highcal"] = "CT: BKG VTX CAL-Med";
    m_TobyToBillMap["gamma/vertex/medcal"]  = "CT: BKG VTX CAL-Med";      
    m_TobyToBillMap["gamma/vertex/thin"]    = "CT: BKG VTX CAL-Low Thin  ";
    m_TobyToBillMap["gamma/vertex/thick"]   = "CT: BKG VTX CAL-Low Thick";
    m_TobyToBillMap["gamma/track/highcal"]  = "CT: BKG 1Tkr CAL-Hi CT ";      
    m_TobyToBillMap["gamma/track/medcal"]   = "CT: BKG 1Tkr CAL-Med";         
    m_TobyToBillMap["gamma/track/thin"]     = "CT: BKG 1Tkr CAL-Low Thin";    
    m_TobyToBillMap["gamma/track/thick"]    = "CT: BKG Bkg 1Tkr CAL-Low Thick";

    // Mapping between the CT output value names in Toby's scheme to IM's
    m_ImToTobyOutMap.clear();
    m_ImToTobyOutMap["Pr(GoodEnergy)"] = "CTgoodCal";
    m_ImToTobyOutMap["Pr(CORE)"]       = "CTgoodPsf";
    m_ImToTobyOutMap["Pr(VTX)"]        = "CTvertex";
    m_ImToTobyOutMap["Pr(GAM)"]        = "CTgamma";
}

const ITreeFactory::ITree& xmlTreeFactory::operator()(const std::string& name)
{
    std::string newName = m_TobyToBillMap[name];

    if (newName != "")
    {
        std::string predEng = "PredictEngineNode";

        // Retrieve vector of PredictEngineNodes 
        std::vector<IImActivityNode*> nodeVec = m_imSheet->getActivityINodeVec(predEng);

        for (std::vector<IImActivityNode*>::iterator nodeIter = nodeVec.begin(); nodeIter != nodeVec.end(); nodeIter++)
        {
            IImActivityNode* iNode = *nodeIter;

            if (iNode->getName() != newName) continue;

            // Convert to a PredictEngineNode
            PredictEngineNode* node = dynamic_cast<PredictEngineNode*>(iNode);

            // Parse the desired tree from the xml document
            DecisionTree* tree = node->getDecisionTree();

            // Retrieve the list of variables used by this tree
            const ImSheetBuilder::StringList& tmpVarNames = node->getInputVarList();

            ImSheetBuilder::StringList& varNames = const_cast<ImSheetBuilder::StringList&>(tmpVarNames);

            GleamValues* gleamVals = new GleamValues(varNames, m_localVals, m_ImToTobyOutMap, m_lookup);

            m_trees.push_back(new Tree(tree, gleamVals));

            break;
        }
    }

    return *m_trees.back();
}


double xmlTreeFactory::Tree::operator()()const
{
    return (*m_dt)(*m_vals);
}
std::string xmlTreeFactory::Tree::title()const
{
    return m_dt->title();
}
xmlTreeFactory::Tree::~Tree()
{
    delete m_dt;
    delete m_vals;
}
xmlTreeFactory::~xmlTreeFactory()
{
    ///@TODO: delete trees
}

xmlTreeFactory::LocalTupleValues::LocalTupleValues(ITreeFactory::ILookupData& lookup) 
{
    m_valsMap.clear();

    // Create the entries in the local variables map
    m_valsMap["BestEnergy"]       = 0.;
    m_valsMap["BestLogEnergy"]    = 0.;
    m_valsMap["BestEnergyProb"]   = 0.;
    m_valsMap["TkrEnergyFrac"]    = 0.;
    m_valsMap["TkrTotalHitsNorm"] = 0.;
    m_valsMap["TkrTotalHitRatio"] = 0.;
    m_valsMap["Tkr12DiffHits"]    = 0.;
    m_valsMap["TkrLATEdge"]       = 0.;
    m_valsMap["PSFEneProbPrd"]    = 0.;
    m_valsMap["CalMaxXtalRatio"]  = 0.;

    // Extract the locations for the ntuple variables
    std::pair<bool, const void*> EvtEnergyCorrPair  = lookup("EvtEnergyCorr");
    std::pair<bool, const void*> TkrEnergyCorrPair  = lookup("TkrEnergyCorr");
    std::pair<bool, const void*> TkrTotalHitsPair   = lookup("TkrTotalHits");
    std::pair<bool, const void*> Tkr1FirstLayerPair = lookup("Tkr1FirstLayer");
    std::pair<bool, const void*> Tkr1HitsPair       = lookup("Tkr1Hits");
    std::pair<bool, const void*> Tkr2HitsPair       = lookup("Tkr2Hits");
    std::pair<bool, const void*> Tkr1X0Pair         = lookup("Tkr1X0");
    std::pair<bool, const void*> Tkr1Y0Pair         = lookup("Tkr1Y0");
    std::pair<bool, const void*> CTgoodPsfPair      = lookup("CTgoodPsf");
    std::pair<bool, const void*> CalXtalMaxEnePair  = lookup("CalXtalMaxEne");
    std::pair<bool, const void*> CalEnergyRawPair   = lookup("CalEnergyRaw");
    
    m_tupleVals["EvtEnergyCorr"]  = (float*)EvtEnergyCorrPair.second;
    m_tupleVals["TkrEnergyCorr"]  = (float*)TkrEnergyCorrPair.second;
    m_tupleVals["TkrTotalHits"]   = (float*)TkrTotalHitsPair.second;
    m_tupleVals["Tkr1FirstLayer"] = (float*)Tkr1FirstLayerPair.second;
    m_tupleVals["Tkr1Hits"]       = (float*)Tkr1HitsPair.second;
    m_tupleVals["Tkr2Hits"]       = (float*)Tkr2HitsPair.second;
    m_tupleVals["Tkr1X0"]         = (float*)Tkr1X0Pair.second;
    m_tupleVals["Tkr1Y0"]         = (float*)Tkr1Y0Pair.second;
    m_tupleVals["CTgoodPsf"]      = (float*)CTgoodPsfPair.second;
    m_tupleVals["CalXtalMaxEne"]  = (float*)CalXtalMaxEnePair.second;
    m_tupleVals["CalEnergyRaw"]   = (float*)CalEnergyRawPair.second;

}

double xmlTreeFactory::LocalTupleValues::getValue(const std::string& name) const
{
    double value = 0.;

    // I would prefer not to do it his way but this is expedient for now...
    //BestEnergy
    if (name == "BestEnergy")
    {
        value = *(m_tupleVals.find("EvtEnergyCorr")->second);

        // Here is what we are supposed to do... currently can't?
        //BestEnergy = ifelse(ProfileProb > ParamProb & ProfileProb > LastLayerProb & 
        //                    ProfileProb > TrackerProb,  CalCfpEnergy, BestEnergy)
    }
    //BestLogEnergy
    else if (name == "BestLogEnergy")
    {
        std::string sBestEnergy = "BestEnergy";
        double BestEnergy = getValue(sBestEnergy);
        value = log(BestEnergy) / 2.3026;
    }
    //BestEnergyProb
    else if (name == "BestEnergyProb")
    {
        value = *(m_tupleVals.find("EvtEnergyCorr")->second);
    }

    //TkrEnergyFrac 
    else if (name == "TkrEnergyFrac")
    {
        value = *(m_tupleVals.find("TkrEnergyCorr")->second) / *(m_tupleVals.find("EvtEnergyCorr")->second);
    }

    //TkrTotalHitsNorm
    else if (name == "TkrTotalHitsNorm")
    {
        value = *(m_tupleVals.find("TkrTotalHits")->second) / sqrt(*(m_tupleVals.find("Tkr1FirstLayer")->second) - 1);
    }

    //TkrTotalHitRatio
    else if (name == "TkrTotalHitRatio")
    {
        value = (*(m_tupleVals.find("Tkr1Hits")->second) + *(m_tupleVals.find("Tkr2Hits")->second))/ *(m_tupleVals.find("TkrTotalHits")->second);
    }

    //Tkr12DiffHits
    else if (name == "Tkr12DiffHits")
    {
        value = *(m_tupleVals.find("Tkr1Hits")->second) - *(m_tupleVals.find("Tkr2Hits")->second);
    }

    //TkrLATEdge
    else if (name == "TkrLATEdge")
    {
        value = 740. - std::max(fabs(*(m_tupleVals.find("Tkr1X0")->second)), fabs(*(m_tupleVals.find("Tkr1Y0")->second)));
    }

    //PSFEneProbPrd
    else if (name == "PSFEneProbPrd")
    {
        double      CTgoodPsf       = *(m_tupleVals.find("CTgoodPsf")->second);
        std::string sBestEnergyProb = "BestEnergyProb";
        double      BestEnergyProb  = getValue(sBestEnergyProb);
        value = sqrt(BestEnergyProb*BestEnergyProb + CTgoodPsf*CTgoodPsf);
    }

    //CalMaxXtalRatio
    else if (name == "CalMaxXtalRatio")
    {
        value = *(m_tupleVals.find("CalXtalMaxEne")->second) / *(m_tupleVals.find("CalEnergyRaw")->second);
    }

    return value;
}

