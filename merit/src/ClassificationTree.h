/** @file ClassificationTree.h
@brief 


$Header$
*/
#ifndef CLASSIFICATIONTREE_H
#define  CLASSIFICATIONTREE_H

#include <string>
#include <vector>
#include "classification/Tree.h"

class Tuple;

class ClassificationTree 
{
public:
    /** set up the tree:
    * @param t The input tuple -- will create a new column with the output
    * @param xml_file 
    */
    ClassificationTree( Tuple&t,  std::string xml_file="");
    /** run the classification
    */
    void execute();  

    ~ClassificationTree();

private:
    classification::Tree * m_classifier;
};

#endif