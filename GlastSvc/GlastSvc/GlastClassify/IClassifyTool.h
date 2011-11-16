/**
 * @class IClassifyTool
 *
 * @brief Interface to the track fitting tools. 
 * Basically, provides an interface to various single track fit 
 * classes. Currently there is but one, but this allows for future
 * ideas/expansion/etc.
 *
 * @author Tracy Usher
 */

#ifndef IClassifyTool_H
#define IClassifyTool_H

#include "GaudiKernel/IAlgTool.h"
#include "GaudiKernel/IProperty.h"

static const InterfaceID IID_IClassifyTool("IClassifyTool", 1 , 0);

namespace GlastClassify
{
    class Item;
};

class IClassifyTool : virtual public IAlgTool
{
public:
    
    /// Typedef the map between variable names and locations
    typedef std::map<std::string, GlastClassify::Item*> VarNameToValueMap;

    /// Retrieve interface ID
    static const InterfaceID& interfaceID() { return IID_IClassifyTool; }

    /// @brief Set up the classification tree analysis for this tool
    virtual StatusCode setUpClassification(VarNameToValueMap& varMap,
                                           std::string&       xmlFileName,
                                           std::string&       tupleName,
                                           std::string        tupleFileName = "") = 0;

    /// @brief Given a pattern track, perform the track fit
    virtual StatusCode runClassification()=0;

    /// @brief Once classification run this will look up and return the value of a given
    ///        variable. If the variable has been found successfully then it returns true
    virtual bool       getVariable(const std::string& varName, GlastClassify::Item* &varItem) = 0;
};
#endif
