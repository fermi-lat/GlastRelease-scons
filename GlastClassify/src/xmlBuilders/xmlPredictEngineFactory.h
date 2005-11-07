/**@file xmlPredictEngineFactory.h

@brief declaration of class xmlPredictEngineFactory
@author T. Usher

$Header$
*/

#ifndef GlastClassify_xmlPredictEngineFactory_h
#define GlastClassify_xmlPredictEngineFactory_h

#include "xmlFactoryBase.h"
#include "IxmlEngineFactory.h"
#include "../ImActivityNodes/IImActivityNode.h"

class DecisionTree;

/** @class xmlPredictEngineFactory
@brief A factory for accessing decision trees

*/
class xmlPredictEngineFactory : public xmlFactoryBase, public IxmlEngineFactory
{
public:
    typedef std::vector<std::string> StringList;

    /** @brief constructotor sets up for tree production 
    @param path file path to a folder containing tree data
    @param lookup Instance of a class supplied to by user, which is called back 
    to find address of each variable
    */
    xmlPredictEngineFactory(std::ostream& log=std::cout, int iVerbosity=0);

    /** @param name a folder name completing the path to the folder containing the tree data   
     @return a reference to a new tree. See also the evaluate() method.
     */
    virtual IImActivityNode* operator()(const DOMElement* element);

    virtual ~xmlPredictEngineFactory();

private:

    /** Build a DecisionTree straight from the top DOMElement
        @param xmlActivityNode is a pointer to the top DOMElement in the forest
    */
    DecisionTree* parseForest(const DOMElement* xmlActivityNode);

    // Mapping between the independent variable in this CT and its "index"
    typedef std::map<std::string, int> tupleVarIndexMap;

    void parseTree(DOMElement* xmlTreeModel, DecisionTree* decisionTree, std::map<std::string, int>& varIndexMap);

    void parseNode(DOMElement* xmlElement, int nodeId, DecisionTree* decisionTree, std::map<std::string, int>& varIndexMap);

    /// Find output variable name for given Classification Tree
    std::string getCTOutputName(const DOMElement* xmlActivityNode);

    /// Find the list of independent variables used by this Classification Tree
    tupleVarIndexMap buildVarIndexMap(DOMElement* xmlNode);

    /// Name of the output variable (e.g. Pr(CORE) or something)
    std::string                  m_outVarName;

    /// "specified category" name - for determining which "yprob" value to use
    std::string                  m_specCatName;

    /// List of variable names used by this tree
    StringList                   m_varNames;

    /// Index to use when extracting value of yprob
    int                          m_yProbIndex;
    int                          m_numProbVals;

    std::ostream&                m_log;         //! output to this stream
    int                          m_outputLevel; //! output level (verbosity)

};


#endif
