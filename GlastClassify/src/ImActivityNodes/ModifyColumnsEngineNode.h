#ifndef ModifyColumnsEngineNode_h
#define ModifyColumnsEngineNode_h

#include "IImActivityNode.h"

#include <vector>
#include <map>
#include <string>
#include <iostream>

/** @class ModifyColumnsEngineNode
*  @brief  Describes an ActivityNode encountered in reading an IM xml file
*
*/

class XTcolumnValBase;

class ModifyColumnsEngineNode : public IImActivityNode
{
public:
    typedef std::map<XTcolumnValBase*,XTcolumnValBase*> XprsnNodeMap;

    ModifyColumnsEngineNode(const std::string& type, const std::string& name, const std::string& id) :
                     m_type(type), m_name(name), m_id(id) {m_nodeMap.clear();}
    virtual ~ModifyColumnsEngineNode() {}

    // 
    virtual const std::string& getType()           const {return m_type;}
    virtual const std::string& getName()           const {return m_name;}
    virtual const std::string& getId()             const {return m_id;}
    virtual const IImActivityNodeMap& getNodeMap() const {return m_nodeMap;}
    
    // Execute the node and its daughters
    virtual void execute();

    virtual void setNodeLink(int port, IImActivityNode* linkToNode) {m_nodeMap[port] = linkToNode;}

    void setXprsnNodeMap(XprsnNodeMap& nodeMap)   {m_xprsnNodeMap = nodeMap;}

    virtual void print(std::ostream& out=std::cout, int depth=0) const;

private:
    std::string        m_type;
    std::string        m_name;
    std::string        m_id;
    IImActivityNodeMap m_nodeMap;

    XprsnNodeMap       m_xprsnNodeMap;
};

#endif // ifdef ModifyColumnsEngineNode_h
