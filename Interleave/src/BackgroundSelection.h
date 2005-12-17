/**


*/


#ifndef BackgroundSelection_h
#define BackgroundSelection_h

#include <string>
class TTree;
class TFile;

/** @class BackgroundSelection
    @brief manage the selection of background events to merge with signal events

*/
class BackgroundSelection
{
public:

    /** @brief ctor
        @param rootFileName name of root file (or files) to open
        @param treename name of the TTree containing the tuple data
    */
    BackgroundSelection(const std::string& rootFileName, const std::string& treename);

    ~BackgroundSelection();

    /** @param maglat the current magnetic latitude
        @return a TTree, with the event select
    */
    TTree* selectEvent(double maglat);

    /**@brief copy the contents of the current tree to another, except for inactive branches 
       @param other the other tree

    */
    void copyTreeData(TTree* other);

    /**@brief the trigger rate for the given magnetic latitude
    */
    double triggerRate(double maglat);

private:
    unsigned int m_event;
    TTree* m_tree;
    TFile* m_file;
};


#endif