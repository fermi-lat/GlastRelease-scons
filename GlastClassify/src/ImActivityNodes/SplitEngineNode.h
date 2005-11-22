#ifndef SplitEngineNode_h
#define SplitEngineNode_h

#include "IImActivityNode.h"

#include <vector>
#include <map>
#include <string>
#include <iostream>

class IXTExprsnNode;

/** @class SplitEngineNode
*  @brief  Describes an ActivityNode encountered in reading an IM xml file
*
*/

class SplitEngineNode : public IImActivityNode
{
public:
    typedef std::vector<std::string> StringList;

    SplitEngineNode(const std::string& type, const std::string& name, const std::string& id) :
                     m_type(type), m_name(name), m_id(id), m_expression("") 
                     {m_nodeMap.clear();}
    virtual ~SplitEngineNode() {}

    // 
    virtual const std::string& getType()           const {return m_type;}
    virtual const std::string& getName()           const {return m_name;}
    virtual const std::string& getId()             const {return m_id;}
    virtual const IImActivityNodeMap& getNodeMap() const {return m_nodeMap;}
    
    // Execute the node and its daughters
    virtual void execute();

    virtual void setNodeLink(int port, IImActivityNode* linkToNode) {m_nodeMap[port] = linkToNode;}

    void setExpression(std::string exp) {m_expression = exp;}
    void setXTExprsnNode(IXTExprsnNode* node) {m_xprsnNode = node;}

    virtual void print(std::ostream& out=std::cout, int depth=0) const;

private:
    std::string        m_type;
    std::string        m_name;
    std::string        m_id;
    IImActivityNodeMap m_nodeMap;

    std::string        m_expression;

    IXTExprsnNode*     m_xprsnNode;
};

#endif // ifdef SplitEngineNode_h
