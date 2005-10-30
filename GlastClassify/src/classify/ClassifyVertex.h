/**@file ClassifyVertex.h
@brief 
$Header$


*/
#ifndef GlastClassify_ClassifyVertex_h
#define GlastClassify_ClassifyVertex_h

#include "classifier/Classifier.h"
#include "GlastClassify.h"

#include <cmath>

class ClassifyVertex : public GlastClassify
{
public:

    ClassifyVertex(const std::string& info_path)
        : GlastClassify ( info_path)
        , McDirErr      ( "McDirErr")
        , McTkr1DirErr  ( "McTkr1DirErr")
    {
    }

    //function to generate good test

    virtual bool isgood()
    {
        return  McDirErr < McTkr1DirErr;
    }

private:
    Entry McDirErr;
    Entry McTkr1DirErr;
};

#endif
