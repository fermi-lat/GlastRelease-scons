#ifndef TEST_ROOTINTERFACE_CXX 
#define TEST_ROOTINTERFACE_CXX 1

#include "rootUtil/RootInterface.h"
#include <string>
#include <iostream>

#include "TH1F.h"
#include "facilities/Util.h"

int main(int argn, char** argc) {
    std::string filename;

    if( argn>1 ) filename = argc[1];

    if (facilities::Util::expandEnvVar(&filename) < 0) {
        std::cout << "Failed to expand env variable in filename" << std::endl;
        return -1;
    }

    // Setup RootInterface
    RootInterface *ri = new RootInterface(filename, "READ");
    std::cout << "You can disreguard these warnings from ROOT" << std::endl;

    // Load ROOT's histogram library - since we know we will need it in this case
    if (ri->loadDynamicLibrary("libHist") < 0) {
        std::cout << "Failed to load libHist" << std::endl;
        return -1;
    }

    // Attempt to retrieve the first histogram located on the top-level directory
    TObject *myObj = ri->retrieveObjByName("myHist");
    if (myObj == 0) {
        std::cout << "Failed to find the first histogram" << std::endl;
        return -1;
    }
    
    std::cout << "Retrieved the first Histogram" << std::endl;

    TH1F *hist = dynamic_cast<TH1F*>(myObj);
    std::string type = myObj->ClassName();

    // Try the retrieveObjByType method for the same histogram
    TObject *myObj2 = ri->retrieveObjByType(type);
    if (myObj2 == 0) {
        std::cout << "Failed to retrieve the first histogram using retrieveObjByType" << std::endl;
        return -1;
    }
    
    // Now attempt to retrieve the second histogram located in a subdirectory
    TObject *myObj3 = ri->retrieveObjByName("myHist2");
    if (myObj3 == 0) {
        std::cout << "Failed to find myHist2" << std::endl;
        return -1;
    }

    std::cout << "Retrieved the second histogram" << std::endl;

    // cleanup
    delete ri;
    ri = 0;

    std::cout << "Completed the test successfully" << std::endl;
    return 0;
}

#endif