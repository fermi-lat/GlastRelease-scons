/** 
 * @file testCheckSum.cxx
 *
 * @brief simple program to check the checkSum class.
 *
 * You need to copy a RootTupleSvc.root to where you are, and run.
 *
 * @author Michael Kuss
 *
 * $Header$
 */

#include "../checkSum.h"

#include "facilities/Util.h"

#include "TFile.h"
#include "TTree.h"
#include "TLeafD.h"
#include "TLeafF.h"

#include <iostream>
#include <iomanip>

#define FATAL   ( 1 )
#define ERROR   ( FATAL )
#define WARNING ( ERROR )
#define INFO    ( WARNING )
// e.g., to disable DEBUG and VERBOSE:
#define DEBUG   ( INFO && 0 )
#define VERBOSE ( DEBUG )

int main() {
    std::string ntupleFileName = "RootTupleSvc.root";
    std::string checkSumFileName = "testCheckSum.txt";

    TFile f(ntupleFileName.c_str());
    TTree* meritTuple = (TTree*)f.Get("MeritTuple");

    const TObjArray* leafCol = meritTuple->GetListOfLeaves();
    const int size = leafCol->GetEntries();
    if ( DEBUG )
        std::cout << "TTree " << meritTuple->GetName() << " has " << size << " leaves "
                  << std::endl;

    checkSum mySum(checkSumFileName);

    for ( int event=0; event<10; ++event ) {
        TTree* newTree = new TTree("newTree", "newTree");

        for ( int item=0; item<size; ++item ) {
            // there exists a TTreeFriendLeafIter, but how to use it?

            TLeaf* leaf = dynamic_cast<TLeaf*>(leafCol->At(item));
            const std::string className = leaf->ClassName();
            const std::string name = leaf->GetName();
            if ( DEBUG )
                std::cout << item << ' ' << name << ' ' << className;

            if ( className == "TLeafD" ) {
                Double_t* value = new Double_t;
                leaf->SetAddress(value);
                meritTuple->GetEvent(event);
                newTree->Branch(name.c_str(), value, (name+"/D").c_str());
                if ( DEBUG )
                    std::cout << ' ' << std::setprecision(25) << *value
                              << std::setprecision(0) << std::endl;
            }
            else if ( className == "TLeafF" ) {
                Float_t* value = new Float_t;
                leaf->SetAddress(value);
                meritTuple->GetEvent(event);
                newTree->Branch(name.c_str(), value, (name+"/F").c_str());
                if ( DEBUG )
                    std::cout << ' ' << std::setprecision(25) << *value
                              << std::setprecision(0) << std::endl;
            }
            else if ( className == "TLeafI" ) {
                Int_t* value = new Int_t;
                leaf->SetAddress(value);
                meritTuple->GetEvent(event);
                newTree->Branch(name.c_str(), value, (name+"/I").c_str());
                if ( DEBUG )
                    std::cout << ' ' << std::setprecision(25) << *value
                              << std::setprecision(0) << std::endl;
            }
            else {
                if ( DEBUG )
                    std::cout << std::endl;
                if ( WARNING )
                    std::cout << "class " << className << " is not implemented!"<<std::endl;
            }
        }

        mySum.write(newTree);
        delete newTree;
    }

    return 0;
}
