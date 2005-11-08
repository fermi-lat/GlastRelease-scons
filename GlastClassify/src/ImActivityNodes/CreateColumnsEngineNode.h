#ifndef CreateColumnsEngineNode_h
#define CreateColumnsEngineNode_h

#include "IImActivityNode.h"

#include <vector>
#include <map>
#include <string>
#include <iostream>

/** @class CreateColumnsEngineNode
*  @brief  Describes an ActivityNode encountered in reading an IM xml file
*
*/

class CreateColumnsEngineNode : public IImActivityNode
{
public:
    typedef std::vector<std::string> StringList;

    CreateColumnsEngineNode(const std::string& type, const std::string& name, const std::string& id) :
                     m_type(type), m_name(name), m_id(id) 
                     {m_nodeVec.clear(); m_columnNames.clear(); 
                      m_columnExpressions.clear(); m_parsedColExps.clear();}
    ~CreateColumnsEngineNode() {}

    // 
    virtual const std::string& getType()           const {return m_type;}
    virtual const std::string& getName()           const {return m_name;}
    virtual const std::string& getId()             const {return m_id;}
    virtual const IImActivityNodeVec& getNodeVec() const {return m_nodeVec;}

    virtual void setNodeLink(IImActivityNode* linkToNode) {m_nodeVec.push_back(linkToNode);}

    void setColumnNames(StringList& colNames)     {m_columnNames       = colNames;}
    void setColumnExpressions(StringList& colExp) {m_columnExpressions = colExp;}
    void setParsedColExps(std::vector<StringList>& colExps) {m_parsedColExps = colExps;}

    virtual void print(std::ostream& out=std::cout, int depth=0) const;

private:
    std::string        m_type;
    std::string        m_name;
    std::string        m_id;
    IImActivityNodeVec m_nodeVec;

    StringList         m_columnNames;
    StringList         m_columnExpressions;
    std::vector<StringList> m_parsedColExps;
};

#endif // ifdef CreateColumnsEngineNode_h
