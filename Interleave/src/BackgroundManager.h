#ifndef BackgroundManager_h
#define BackgroundManager_h 

#include <string>
#include <vector>
#include <map>

class IBkgndTupleSelectTool;
/** 
* @class BackgroundManager
*
* @brief Singleton object providing management and utility routines for
*        background sources used by the interleave. Specifically, it is 
*        designed to enable communication between the SampledBackground
*        class (which is what defines a "source" to FluxSvc) and 
*        BackgroundSelection (the class which takes care of the actual 
*        input of an interleave event). 
*
* November 7, 2007 Tracy Usher
*
* @author Tracy Usher
*
*/

//Class definition
class BackgroundManager
{
public: 
    typedef std::vector<std::string>                      BkgndSrcNameVec;
    typedef std::map<std::string, IBkgndTupleSelectTool*> BkgndNameToSelMap;

  /// Method to return pointer to the singleton object (and instantiate if not
  /// already)
    static  BackgroundManager* instance();

    /// Methods for interfacing to the class
    ///! access rate according to section variable
    ///! @param tvar Name of key variable: if blank, sum all rates
    double downlinkRate(const std::string& tvar="");
    double triggerRate(const std::string& tvar="");

    ///! during initialization, add this source name
    void defineSource(const std::string& tvar);

    ///! These provide access to vectors/maps defining the sources
    BkgndSrcNameVec&   getBkgndSrcNameVec()    {return m_sourceNames;}
    BkgndNameToSelMap& getBkgndNameToSelMap()  {return m_bkgndToSelMap;}
    BkgndNameToSelMap& getSrcNameToSelMap()    {return m_srcToSelMap;}
  
private:
    /// Private methods
    BackgroundManager();
   ~BackgroundManager();

    /// Private data
    /// Pointer to the class to make it a singleton
    static BackgroundManager* m_instance;

    double                    m_defaultRate;

    BkgndSrcNameVec           m_sourceNames;
    BkgndNameToSelMap         m_bkgndToSelMap;
    BkgndNameToSelMap         m_srcToSelMap;
};

#endif




