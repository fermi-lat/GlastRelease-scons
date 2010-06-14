/**@file main.cxx
@brief main program for application to test GLAST classification trees

*/

#include "GlastClassify/TreeFactory.h"
#include <stdexcept>
#include <iostream>
#include <fstream>

#include <string>
#include <vector>
using namespace GlastClassify;

static std::vector<double> values;
class TestLookup : public TreeFactory::ILookupData {
public:
    TestLookup(){}
    const double * operator()(const std::string& name){
        std::cout << "Looking up: " << name << std::endl;
        values.push_back(values.size());
        return &values.back();
    }
};


class ClassifyTest {
public:
    ClassifyTest( const std::string& info_path, const std::string& name, std::ostream& inlog)
        : m_log(&inlog)
    {
        log() << "Creating the factory:" << std::endl;

        TreeFactory factory(info_path, TestLookup());

        log() << "Ask factory for a tree: " << std::endl;
        const TreeFactory::Tree& goodcal= factory(name);
        log() << "Tree's title: " << goodcal.title() << std::endl;

        log() << name  << " value from the tree: " << goodcal() << std::endl;
    }
  
private:
    std::ostream * m_log;
    std::ostream& log(){return * m_log;}
};



int main(int argc, char ** argv)
{
    int rc=0;
    std::string name("goodcal");
    try {
        std::cout << "Gleam classification test:\n\ttree:\t"
            << name << std::endl;
        std::string tree_path(std::string(::getenv("GLASTCLASSIFYROOT"))+"/data");

        std::cout << "\tpath:\t" << tree_path << std::endl;

        ClassifyTest( tree_path, name, std::cout);

    }catch (const std::exception & error)    {
        std::cerr << "Caught exception \"" << error.what() << "\"" <<std::endl;
        rc= 1;
    }
    return rc;
}
