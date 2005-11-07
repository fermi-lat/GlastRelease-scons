#ifndef SplitEngineNode_h
#define SplitEngineNode_h

#include "IImActivityNode.h"

#include <vector>
#include <map>
#include <string>
#include <iostream>

/** @class SplitEngineNode
*  @brief  Describes an ActivityNode encountered in reading an IM xml file
*
*/

class SplitEngineNode : public IImActivityNode
{
public:

    SplitEngineNode(const std::string& type, const std::string& name, const std::string& id) :
                     m_type(type), m_name(name), m_id(id), m_expression("") {m_nodeVec.clear();}
    ~SplitEngineNode() {}

    // 
    virtual const std::string& getType()           const {return m_type;}
    virtual const std::string& getName()           const {return m_name;}
    virtual const std::string& getId()             const {return m_id;}
    virtual const IImActivityNodeVec& getNodeVec() const {return m_nodeVec;}

    virtual void setNodeLink(IImActivityNode* linkToNode) {m_nodeVec.push_back(linkToNode);}

    void setExpression(std::string exp) {m_expression = exp;}

    virtual void print(std::ostream& out=std::cout, int depth=0) const;

private:
    std::string        m_type;
    std::string        m_name;
    std::string        m_id;
    IImActivityNodeVec m_nodeVec;

    std::string        m_expression;
};

#endif // ifdef SplitEngineNode_h
