#ifndef ReadTextFileEngineNode_h
#define ReadTextFileEngineNode_h

#include "IImActivityNode.h"

#include <vector>
#include <map>
#include <string>
#include <iostream>

/** @class ReadTextFileEngineNode
*  @brief  Describes an ActivityNode encountered in reading an IM xml file
*
*/
template <class T> class XTcolumnVal;

class ReadTextFileEngineNode : public IImActivityNode
{
public:
    typedef std::vector<std::string> StringList;

    ReadTextFileEngineNode(const std::string& type, const std::string& name, const std::string& id) :
                     m_type(type), m_name(name), m_id(id) {m_nodeMap.clear();}
    virtual ~ReadTextFileEngineNode() {}

    // 
    virtual const std::string& getType()           const {return m_type;}
    virtual const std::string& getName()           const {return m_name;}
    virtual const std::string& getId()             const {return m_id;}
    virtual const IImActivityNodeMap& getNodeMap() const {return m_nodeMap;}
    
    // Execute the node and its daughters
    virtual void execute();

    virtual void setNodeLink(int port, IImActivityNode* linkToNode) {m_nodeMap[port] = linkToNode;}

    virtual void print(std::ostream& out=std::cout, int depth=0) const;

    void setXtColumnVal(XTcolumnVal<double>* xtColumnVal) {m_WriteTupleRow = xtColumnVal;}

private:
    std::string          m_type;
    std::string          m_name;
    std::string          m_id;
    XTcolumnVal<double>* m_WriteTupleRow;
    IImActivityNodeMap   m_nodeMap;
};

#endif // ifdef ReadTextFileEngineNode_h
