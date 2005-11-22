/**@file xmlCreateColumnsEngineFactory.h

@brief declaration of class xmlCreateColumnsEngineFactory
@author T. Usher

$Header$
*/

#ifndef GlastClassify_xmlCreateColumnsEngineFactory_h
#define GlastClassify_xmlCreateColumnsEngineFactory_h

#include "xmlFactoryBase.h"
#include "IxmlEngineFactory.h"
#include "src/ImActivityNodes/IImActivityNode.h"

class DecisionTree;

/** @class xmlCreateColumnsEngineFactory
@brief A factory for accessing decision trees

*/
class xmlCreateColumnsEngineFactory : public xmlFactoryBase , public IxmlEngineFactory
{
public:
    /** @brief constructotor sets up for tree production 
    @param path file path to a folder containing tree data
    @param lookup Instance of a class supplied to by user, which is called back 
    to find address of each variable
    */
    xmlCreateColumnsEngineFactory(XTExprsnParser& parser);

    /** @param name a folder name completing the path to the folder containing the tree data   
     @return a reference to a new tree. See also the evaluate() method.
     */
    IImActivityNode* operator()(const DOMElement* element);

    virtual ~xmlCreateColumnsEngineFactory();

private:
};


#endif
