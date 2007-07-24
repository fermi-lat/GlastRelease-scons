/**@file xmlFindOutputVars.h

@brief declaration of class xmlFindOutputVars
@author T. Usher

$Header$
*/

#ifndef GlastClassify_xmlFindOutputVars_h
#define GlastClassify_xmlFindOutputVars_h

#include "xmlFactoryBase.h"
#include "IxmlEngineFactory.h"
#include "src/ImActivityNodes/IImActivityNode.h"

/** @class xmlFindOutputVars
@brief A factory for accessing decision trees

*/
class xmlFindOutputVars : public xmlFactoryBase
{
public:
    /** @brief constructotor sets up for tree production 
    @param path file path to a folder containing tree data
    @param lookup Instance of a class supplied to by user, which is called back 
    to find address of each variable
    */
    xmlFindOutputVars(XTExprsnParser& parser);

    /** @param name a folder name completing the path to the folder containing the tree data   
     @return a reference to a new tree. See also the evaluate() method.
     */
    int operator()(const DOMElement* element);

    virtual ~xmlFindOutputVars();

private:
    /** @brief finds variables defined in createColumnsEngineNodes
    */
    int numCreateColumnsVars(const DOMElement* element);

    /** @brief finds variables output from predictEngineNodes
    */
    int numPredictEngineVars(const DOMElement* element);
};


#endif
