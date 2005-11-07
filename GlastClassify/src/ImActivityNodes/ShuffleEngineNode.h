#ifndef ShuffleEngineNode_h
#define ShuffleEngineNode_h

#include "IImActivityNode.h"

#include <vector>
#include <map>
#include <string>
#include <iostream>

/** @class ShuffleEngineNode
*  @brief  Describes an ActivityNode encountered in reading an IM xml file
*
*/

class ShuffleEngineNode : public IImActivityNode
{
public:
    typedef std::vector<std::string> StringList;

    ShuffleEngineNode(const std::string& type, const std::string& name, const std::string& id) :
                     m_type(type), m_name(name), m_id(id) {m_nodeVec.clear();}
    ~ShuffleEngineNode() {}

    // 
    virtual const std::string& getType()           const {return m_type;}
    virtual const std::string& getName()           const {return m_name;}
    virtual const std::string& getId()             const {return m_id;}
    virtual const IImActivityNodeVec& getNodeVec() const {return m_nodeVec;}

    virtual void setNodeLink(IImActivityNode* linkToNode) {m_nodeVec.push_back(linkToNode);}

    virtual void print(std::ostream& out=std::cout, int depth=0) const;

private:
    std::string        m_type;
    std::string        m_name;
    std::string        m_id;
    IImActivityNodeVec m_nodeVec;
};

#endif // ifdef ShuffleEngineNode_h
