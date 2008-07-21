/**@file main.cxx
@brief main program for application to test GLAST classification trees

$Header$

*/

#include "GlastClassify/TreeFactory.h"
#include "facilities/commonUtilities.h"
#include <stdexcept>
#include <iostream>

#include <string>
#include <vector>

using namespace GlastClassify;

static std::vector<float> values;

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// example of a ITem subclass that dereferences a pointer
class TestItem : public Item {
public:
    TestItem(float * val): m_value(val){}
    operator double()const { return *m_value;}
private:
    float * m_value;
};
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class TestTuple : public ITupleInterface{
public:
    const Item * getItem(const std::string& name)const
    {
        std::cout << "Looking up: " << name << std::endl;
        float value = values.size(); // very klugy, just exercise logic
        values.push_back( value);
        const TestItem* item = new TestItem(&values.back());
        return item;
    }
    void addItem(const std::string&, float & ){}
    void addItem(const std::string&, double & ){}
};
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


class ClassifyTest {
public:
    ClassifyTest( const std::string& info_path, const std::string& name, std::ostream& inlog)
        : m_log(&inlog)
    {
        log() << "Creating the factory:" << std::endl;

        ITupleInterface& tuple = *new TestTuple();
        TreeFactory::LocalDictionary dict; // empty
        TreeFactory factory(info_path, tuple, dict);

        log() << "Ask factory for a tree: " << std::endl;
        const TreeFactory::Tree& goodcal= factory(name);
        log() << "Tree's title: " << goodcal.title() << std::endl;
        

        log() << name  << " value from the tree: " << goodcal() << std::endl;
    }
  
private:
    std::ostream * m_log;
    std::ostream& log(){return * m_log;}
};

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

int main(int, char ** )
{
    facilities::commonUtilities::setupEnvironment();
    int rc=0;
    std::string name("goodcal");
    try {
        std::cout << "Gleam classification test:\n\ttree:\t"
            << name << std::endl;
        const char * root= ::getenv("GLASTCLASSIFYROOT");
	std::string tree_path = facilities::commonUtilities::getDataPath("GlastClassify");

        std::cout << "\tpath:\t" << tree_path << std::endl;

        ClassifyTest( tree_path, name, std::cout);

    }catch (const std::exception & error)    {
        std::cerr << "Caught exception \"" << error.what() << "\"" <<std::endl;
        rc= 1;
    }
    return rc;
}
