/** @file BackgroundSelection.h

    @brief declaration of the BackGroundSelection class

$Header$

*/


#ifndef BackgroundSelection_h
#define BackgroundSelection_h

#include <string>
#include <vector>
class TTree;
class TFile;
class TLeaf;

/** @class BackgroundSelection
    @brief manage the selection of background events to merge with signal events
    @author Dan Flath

It uses the name of a tuple variable as a key to look up trigger and downlink rates of a background
source depending on the value of that variable. See the corresponding member functions. Given such
a value, it expects to find a source of such events generated with that value which can be inserted 
into the Gleam output tuple.

*/
class BackgroundSelection
{
public:

    /** @brief ctor
        @param varname tuple variable name
        @param rootFileDirectory name of the directory containing the root files, or a root file for testing
        @param disableList list of strings specifying branches to disable
        @param outputTree pointer to tree to copy to
    */
    BackgroundSelection(const std::string& varname, const std::string& rootFileDirectory,
        std::vector<std::string> disableList,
        TTree* outputTree);

    ~BackgroundSelection();

    ///! The current value of the quantity that we are selecting on
    double value()const;

    /** @brief select an event and copy the contents to the output tree
    */
    void selectEvent();

 
    /**@brief the current trigger rate for the current value
    */
    double triggerRate()const;

    /**@brief the downlink rate for the current value
    */
    double downlinkRate()const;

private:
    /** @brief disable branches in the tree, from the disableList
    */
    void disableBranches(TTree*);

    /**@brief Set the addresses of active leaves, so that a GetEvent will perform a copy
    */
    void setLeafPointers(TTree*);

    /**@brief Retrieve the correct tree for current geomagnetic latitude so copyEvent comes from correct sample bin
    */
    void setCurrentTree();

    std::string m_varname; ///< name of the variable that we use to key the rates
    TLeaf* m_varleaf;      ///< corresponding TLeaf for access to current value

    unsigned int m_event, m_eventOffset;

    bool m_singleFileMode; 
    
    TTree* m_outputTree, ///< the tree to copy data to
      *m_inputTree; ///< the tree to copy data from
    
    TFile* m_inputFile; ///< the file that contains the input tree
    std::vector<std::string>m_disableList; ///< list of strings for disable
    std::string m_rootFileDirectory;  ///< either path to a ROOT file to sample from, or a directory containing files

};


#endif
