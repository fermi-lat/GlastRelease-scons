#ifndef newPredictEngineNode_h
#define newPredictEngineNode_h

#include "IImActivityNode.h"

#include <vector>
#include <map>
#include <string>
#include <iostream>

class IXTExprsnNode;
template <class T> class XTcolumnVal;

/** @class newPredictEngineNode
*  @brief  Describes an ActivityNode encountered in reading an IM xml file
*
*/

class newPredictEngineNode : public IImActivityNode
{
public:
    typedef std::vector<std::string>          StringList;
    typedef std::pair<IXTExprsnNode*, double> TreePair;
    typedef std::vector<TreePair>             TreePairVector;

    newPredictEngineNode(const std::string& type, const std::string& name, const std::string& id);
    virtual ~newPredictEngineNode() {}

    // 
    virtual const std::string& getType()           const {return m_type;}
    virtual const std::string& getName()           const {return m_name;}
    virtual const std::string& getId()             const {return m_id;}
    virtual const IImActivityNodeMap& getNodeMap() const {return m_nodeMap;}
    
    // Execute the node and its daughters
    virtual void execute();

    const StringList&  getOutputVarList()          const {return m_outputVar;}
    const StringList&  getInputVarList()           const {return m_inputVar;}

    virtual void setNodeLink(int port, IImActivityNode* linkToNode) {m_nodeMap[port] = linkToNode;}
    void setInputVar(const StringList& inVars)            {m_inputVar = inVars;}
    void setOutputVar(const StringList& outVars)          {m_outputVar = outVars;}
    void addOutputVar(const std::string& outVar)          {m_outputVar.push_back(outVar);}
    void setXTcolumnVal(XTcolumnVal<double>* colVal)      {m_xtColumnVal = colVal;}
    void setPredictVal(XTcolumnVal<std::string>* predict) {m_predict = predict;}
    void setTreePairVector(TreePairVector& treeVec)       {m_trees = treeVec;}
    void addTree(TreePair& treePair)                      {m_trees.push_back(treePair);}

    virtual void print(std::ostream& out=std::cout, int depth=0) const;

private:
    std::string               m_type;
    std::string               m_name;
    std::string               m_id;
    IImActivityNodeMap        m_nodeMap;
    StringList                m_outputVar;
    StringList                m_inputVar;
    XTcolumnVal<double>*      m_xtColumnVal;
    XTcolumnVal<std::string>* m_predict;
    TreePairVector            m_trees;
};

// Need to share a Classification Tree output class here
class CTOutPut
{
public:
    CTOutPut() : m_id(0), m_score(""), m_recordcount(0), m_group(0), m_deviance(0), m_entropy(0), m_gini(0), m_risk(0), m_yprob(0) {};
    CTOutPut(int id, std::string& score, int recordcount, int group, double deviance, double entropy, double gini, int risk, double yprob) :
            m_id(id),
            m_score(score),
            m_recordcount(recordcount),
            m_group(group),
            m_deviance(deviance),
            m_entropy(entropy),
            m_gini(gini),
            m_risk(risk),
            m_yprob(yprob) {};
    ~CTOutPut() {};

    const std::string& getScore() const {return m_score;}
    double             getYProb() const {return m_yprob;}

    friend std::ostream& operator <<(std::ostream& stream, const CTOutPut& node);

private:
    int         m_id;
    std::string m_score;
    int         m_recordcount;
    int         m_group;
    double      m_deviance;
    double      m_entropy;
    double      m_gini;
    int         m_risk;
    double      m_yprob;
};

#endif // ifdef newPredictEngineNode_h
