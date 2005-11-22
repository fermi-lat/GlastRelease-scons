/**@file xmlNewPredictEngine.h

@brief declaration of class newPredictEngine
@author T. Usher

$Header$
*/

#ifndef GlastClassify_newPredictEngine_h
#define GlastClassify_newPredictEngine_h

#include "xmlFactoryBase.h"
#include "IxmlEngineFactory.h"
#include "../ImActivityNodes/newPredictEngineNode.h"

/** @class newPredictEngine
@brief A factory for accessing decision trees

*/
class xmlNewPredictEngineFactory : public xmlFactoryBase, public IxmlEngineFactory
{
public:
    typedef std::vector<std::string> StringList;

    /** @brief constructotor sets up for tree production 
    @param path file path to a folder containing tree data
    @param lookup Instance of a class supplied to by user, which is called back 
    to find address of each variable
    */
    xmlNewPredictEngineFactory(XTExprsnParser& parser);

    /** @param name a folder name completing the path to the folder containing the tree data   
     @return a reference to a new tree. See also the evaluate() method.
     */
    virtual IImActivityNode* operator()(const DOMElement* element);

    virtual ~xmlNewPredictEngineFactory();

private:

    /** Build a DecisionTree straight from the top DOMElement
        @param xmlActivityNode is a pointer to the top DOMElement in the forest
    */
    newPredictEngineNode::TreePairVector parseForest(const DOMElement* xmlActivityNode);

    IXTExprsnNode* parseTree(DOMElement* xmlTreeModel);

    IXTExprsnNode* parseNode(DOMElement* xmlElement);

    /// Find output variable name for given Classification Tree
    std::string getCTOutputName(const DOMElement* xmlActivityNode);

    /// Name of the output variable (e.g. Pr(CORE) or something)
    std::string                  m_outVarName;

    /// "specified category" name - for determining which "yprob" value to use
    std::string                  m_specCatName;

    /// List of variable names used by this tree
    StringList                   m_varNames;

    /// Index to use when extracting value of yprob
    int                          m_yProbIndex;
    int                          m_numProbVals;
};


#endif
