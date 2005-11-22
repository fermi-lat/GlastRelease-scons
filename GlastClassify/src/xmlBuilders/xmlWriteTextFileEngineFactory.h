/**@file xmlWriteTextFileEngineFactory.h

@brief declaration of class xmlWriteTextFileEngineFactory
@author T. Usher

$Header$
*/

#ifndef GlastClassify_xmlWriteTextFileEngineFactory_h
#define GlastClassify_xmlWriteTextFileEngineFactory_h

#include "xmlFactoryBase.h"
#include "IxmlEngineFactory.h"
#include "../ImActivityNodes/IImActivityNode.h"

class DecisionTree;

/** @class xmlWriteTextFileEngineFactory
@brief A factory for accessing decision trees

*/
class xmlWriteTextFileEngineFactory : public xmlFactoryBase, public IxmlEngineFactory 
{
public:
    typedef std::vector<std::string> StringList;

    /** @brief constructotor sets up for tree production 
    @param path file path to a folder containing tree data
    @param lookup Instance of a class supplied to by user, which is called back 
    to find address of each variable
    */
    xmlWriteTextFileEngineFactory(XTExprsnParser& parser);

    /** @param name a folder name completing the path to the folder containing the tree data   
     @return a reference to a new tree. See also the evaluate() method.
     */
    IImActivityNode* operator()(const DOMElement* element);

    virtual ~xmlWriteTextFileEngineFactory();

private:
};


#endif
