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

/** @class BackgroundSelection
    @brief manage the selection of background events to merge with signal events

*/
class BackgroundSelection
{
public:

    /** @brief ctor
        @param rootFileDirectory name of the directory containing the root files
        @param disableList list of strings specifying branches to disable
        @param outputTree pointer to tree to copy to
    */
    BackgroundSelection(const std::string& rootFileDirectory,
        std::vector<std::string> disableList,
        TTree* outputTree);

    ~BackgroundSelection();


    /** @brief select an event and copy the contents to the output tree
        @param maglat the current magnetic latitude
    */
    void selectEvent(double maglat);

 
    /**@brief the trigger rate for the given magnetic latitude
    */
    static double triggerRate(double maglat);

    /**@brief the downlink rate for the given magnetic latitude
    */
    static double downlinkRate(double maglat);

private:
    /** @brief disable branches in the tree, from the disableList
    */
    void disableBranches(TTree*);

    /**@brief Set the addresses of active leaves, so that a GetEvent will perform a copy
    */
    void setLeafPointers();
    void setLeafPointers(TTree*);

    /**@brief Retrieve the correct tree for current geomagnetic latitude so copyEvent comes from correct sample bin
    */
    void setCurrentTree(double maglat);

    unsigned int m_event, m_treeInitialized, m_eventOffset;
    double m_currentBinLower;
    TTree* m_outputTree, ///< the tree to copy data to
      *m_inputTree; ///< the tree to copy data from
    TFile* m_inputFile; ///< the file that contains the input tree
    std::vector<std::string>m_disableList; ///< list of strings for disable
    std::string m_rootFileDirectory;
};


#endif
