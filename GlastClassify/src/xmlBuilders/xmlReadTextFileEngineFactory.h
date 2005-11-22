/**@file xmlReadTextFileEngineFactory.h

@brief declaration of class xmlReadTextFileEngineFactory
@author T. Usher

$Header$
*/

#ifndef GlastClassify_xmlReadTextFileEngineFactory_h
#define GlastClassify_xmlReadTextFileEngineFactory_h

#include "xmlFactoryBase.h"
#include "IxmlEngineFactory.h"
#include "../ImActivityNodes/IImActivityNode.h"

class DecisionTree;

/** @class xmlReadTextFileEngineFactory
@brief A factory for accessing decision trees

*/
class xmlReadTextFileEngineFactory : public xmlFactoryBase, public IxmlEngineFactory 
{
public:
    typedef std::vector<std::string> StringList;

    /** @brief constructotor sets up for tree production 
    @param path file path to a folder containing tree data
    @param lookup Instance of a class supplied to by user, which is called back 
    to find address of each variable
    */
    xmlReadTextFileEngineFactory(XTExprsnParser& parser);

    /** @param name a folder name completing the path to the folder containing the tree data   
     @return a reference to a new tree. See also the evaluate() method.
     */
    IImActivityNode* operator()(const DOMElement* element);

    virtual ~xmlReadTextFileEngineFactory();

private:
};


#endif
