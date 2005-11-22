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
    void setTreePairVector(TreePairVector& treeVec)       {m_trees = treeVec;}
    void addTree(TreePair& treePair)                      {m_trees.push_back(treePair);}

    virtual void print(std::ostream& out=std::cout, int depth=0) const;

private:
    std::string          m_type;
    std::string          m_name;
    std::string          m_id;
    IImActivityNodeMap   m_nodeMap;
    StringList           m_outputVar;
    StringList           m_inputVar;
    XTcolumnVal<double>* m_xtColumnVal;
    TreePairVector       m_trees;
};

#endif // ifdef newPredictEngineNode_h
