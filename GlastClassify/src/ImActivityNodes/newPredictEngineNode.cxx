/** @file ImSheetBuilder.cxx
 *    @brief implementation of classification::Tree; declaration and implementation or its private helper classification::ImSheetBuilder::Node
 *
 *    $Header$
 */

#include "newPredictEngineNode.h"
#include "src/XT/XTtupleVars.h"
#include "src/XT/XprsnTree.h"

#include <iostream>
#include <iomanip>
#include <algorithm>
#include <sstream>
#include <set> 


namespace 
{
    std::string indent(int depth)
    {
        std::string ret("  ");; 
        for( int i =0; i < depth; ++i) ret += "  ";
        return ret;
    }
} // anonymous namespace

newPredictEngineNode::newPredictEngineNode(const std::string& type, const std::string& name, const std::string& id) : 
                                           m_type(type), 
                                           m_name(name), 
                                           m_id(id)  
{
    m_nodeMap.clear();
    m_inputVar.clear();
    m_outputVar.clear();
    m_trees.clear();
}

// Does the "real" work... 
void newPredictEngineNode::print(std::ostream& out, int depth) const
{
    // Output our node ID, type and name
    out << indent(depth) << "ID: " << m_id << ", Type: " << m_type << ", Label: " << m_name << std::endl;

    // What do we set depth to?
    depth = m_nodeMap.size() > 1 ? depth + 1 : depth;

    // Now follow through with all the nodes we point to
    for(IImActivityNodeMap::const_iterator nodeIter = m_nodeMap.begin(); nodeIter != m_nodeMap.end(); nodeIter++)
    {
        nodeIter->second->print(out, depth);
    }

    return;
}

// Define an object to keep track of tree results
class PredNodeInfo
{
public:
    PredNodeInfo() : m_count(0), m_catVarName(""), m_weightSum(0.), m_probSum(0.) {};
    PredNodeInfo(int count, std::string& catVarName, double weightSum, double probSum) :
                 m_count(count),
                 m_catVarName(catVarName),
                 m_weightSum(weightSum),
                 m_probSum(probSum) {};
    ~PredNodeInfo() {};

    int                getCount()      const {return m_count;}
    const std::string& getCatVarName() const {return m_catVarName;}
    double             getWeightSum()  const {return m_weightSum;}
    double             getProbSum()    const {return m_probSum;}

    void   incCount()                {m_count++;}
    void   addWeight(double weight)  {m_weightSum += weight;}
    void   addProb(double prob)      {m_probSum   += prob;}

    void   setCatVarName(std::string& catVarName) {m_catVarName = catVarName;}

private:
    int         m_count;
    std::string m_catVarName;
    double      m_weightSum;
    double      m_probSum;
};

typedef std::map<std::string,PredNodeInfo>           PredNodeMap;
typedef std::map<std::string,PredNodeInfo>::iterator PredNodeMapItr;

// Does the "real" work... 
void newPredictEngineNode::execute()
{
    // Loop through the forest evaluating the trees
    std::string predict = "";

    PredNodeMap predNodeMap;

    for(TreePairVector::iterator treeItr = m_trees.begin(); treeItr != m_trees.end(); treeItr++)
    {
        IXTExprsnNode* tree = treeItr->first;

        //double result = *(reinterpret_cast<const double*>((*tree)()));
        const CTOutPut* ctOutPut = reinterpret_cast<const CTOutPut*>((*tree)());

        double result = ctOutPut->getYProb();

        predict = ctOutPut->getScore();

        PredNodeInfo& nodeInfo = predNodeMap[predict];

        nodeInfo.incCount();
        nodeInfo.addProb(treeItr->second * result);
        nodeInfo.addWeight(treeItr->second);
        nodeInfo.setCatVarName(predict);

        // For printing during running, leaving here to remember how to do it
        //const XTIfElseNode<double>& ifElseNode = dynamic_cast<const XTIfElseNode<double>& >(*tree);
        //std::cout << "Tree #" << treeNum++ << " evaluates to: " << result << std::endl;
        //ifElseNode.printExp(std::cout);

        //weightSum += treeItr->second;
        //runningWght += treeItr->second * result;
    }

    //std::cout << "Tree " << m_name << ", running sum=" << runningWght 
    //          << ", weightSum=" << weightSum;

    PredNodeMapItr predNodeMapItr = predNodeMap.begin();
    double         runningWght    = 0.;
    double         weightSum      = 0.;

    if (predNodeMap.size() > 0)
    {
        int maxCount = -1;
        for(PredNodeMapItr itr = predNodeMap.begin(); itr != predNodeMap.end(); itr++)
        {
            if (itr->second.getCount() > maxCount)
            {
                predNodeMapItr = itr;
                maxCount       = itr->second.getCount();
            }

            // Accumulate weights
            runningWght += itr->second.getProbSum();
            weightSum   += itr->second.getWeightSum();
        }
    }

    PredNodeInfo& predNode = predNodeMapItr->second;

    if (weightSum > 0.) runningWght /= weightSum;

    if (predNode.getCount() > 0)
    {
        //runningWght = predNode.getProbSum() / predNode.getWeightSum();
        predict     = predNode.getCatVarName();
    }

    //std::cout << " gives: " << runningWght << std::endl;
    
    m_xtColumnVal->setDataValue(runningWght);
    m_predict->setDataValue(predict);
    
    // Now follow through with all the daughter nodes we point to
    for(IImActivityNodeMap::const_iterator nodeIter = m_nodeMap.begin(); nodeIter != m_nodeMap.end(); nodeIter++)
    {
        nodeIter->second->execute();
    }

    return;
}
// Override the ostream << operator. Keeping here as an example....
std::ostream& operator<<(std::ostream& stream, const CTOutPut& node)
{
    stream << node.getScore() << "=" << node.getYProb();
    return stream;
}
