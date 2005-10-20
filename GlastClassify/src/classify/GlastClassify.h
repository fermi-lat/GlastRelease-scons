/**@file GlastClassify.h
@brief 
$Header$

*/
#ifndef GlastClassify_GlastClassify_h
#define GlastClassify_GlastClassify_h


#include "classifier/Classifier.h"
#include "classifier/TrainingInfo.h"


class GlastClassify
{
public:

    typedef enum{ALL, ODD, EVEN} Subset;

    /// ctor
    GlastClassify(const std::string& info_path, bool mixed=true);

    virtual ~GlastClassify(){} 

    /**
    @brief operate on the given subset
    @param max_events limit to number of events
    @param set     subset of each file to use
    @param train   train if true, evaluate if false

    */
    void run( unsigned int max_events=0, Subset set=ALL);

    static void setPaths(std::string rootpath, std::string treepath){
        s_rootpath = rootpath+"/";
        s_treepath = treepath+"/";
    }
    /// path to the root files
    static std::string s_rootpath;
    /// path to tree data, input and output
    static std::string s_treepath;

    static bool s_train;

protected:
    /// subclasses may implement this to define the good, or signal events
    /// default value is controlled by setbkgnd(v), v=true for background, false for signal
    virtual bool isgood() {return !m_nobkgnd;};

    /// subclass may override
    virtual void define(std::vector<std::string>& /*all_names*/){};

    /// acceptance cut applied to events in training sample: subclass may override
    virtual bool accept(){return true;}

    /// return index of the variable in the training list (exception if not found)
    int find_index(const std::string& name);

    /// add an entry to the list for making cuts to define sample to train on
    int add_index(const std::string& name);

    /// access to values in the current event
    float datum(int index)const{return (*m_row)[index];};

    void setbkgnd (bool v=false){ m_nobkgnd=v;}; 

    class Entry{
    public:
        Entry( const std::string& name="");
        operator double()const{ return m_cl->datum(m_index);}
    private:
        GlastClassify* m_cl;
        int m_index;
    };
    
    
private:

    void load( unsigned int max_events, Subset set);

    void load(TrainingInfo::StringList input, unsigned int max_events=0, bool good=true);

    void classify();
    void test();

    void current_time(std::ostream& out=std::cout);

    std::string m_title;
    std::ostream * m_log;
    std::ostream& log(){return * m_log;}
    Classifier::Table m_data;
    TrainingInfo m_info;
    const std::vector<float>* m_row; ///< current row in interating thru tuple
    bool m_nobkgnd;
    std::vector<std::string> m_all_names;
    bool m_mixed;       ///< true in file has mixed good/bad: uses isGood() in this case

};

#endif
