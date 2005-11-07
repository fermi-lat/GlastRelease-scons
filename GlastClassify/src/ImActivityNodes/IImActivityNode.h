#ifndef IImActivityNode_h
#define IImActivityNode_h

/** @class IImActivityNode
*  @brief  Abstract interface to define an "Activity Node" in an analysis
*
*/

#include <vector>
#include <string>
#include <iostream>

class IImActivityNode
{
public:
    // Typedef for the a vector of activity nodes
    typedef std::vector<IImActivityNode*> IImActivityNodeVec;

    // Minimal set of virtual functions
    virtual const std::string& getType()           const = 0;  // Node type
    virtual const std::string& getName()           const = 0;  // Node name
    virtual const std::string& getId()             const = 0;  // Node ID
    virtual const IImActivityNodeVec& getNodeVec() const = 0;  // Vector of child nodes
    
    virtual void setNodeLink(IImActivityNode* linkToNode) = 0;

    // output
    virtual void print(std::ostream& out=std::cout, int depth=0) const = 0;
};

#endif // ifdef IImActivityNode_h
