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
        , m_isThin      ( info_path.find("thin")>0)
        , VtxAngle      ( "VtxAngle")
        , Tkr1FirstLayer( "Tkr1FirstLayer")
        , McDirErr      ( "McDirErr")
        , McTkr1DirErr  ( "McTkr1DirErr")
        , CTgoodCal     ( "CTgoodCal")
    {
    }

    //function to generate good test

    virtual bool isgood()
    {
        return  McDirErr < McTkr1DirErr;
    }

    virtual bool accept()
    {
        if( CTgoodCal< 0.25 || VtxAngle==0) return false;

        bool thin(Tkr1FirstLayer > 5);
        return thin==m_isThin;
    }

private:
    Entry VtxAngle;      
    Entry Tkr1FirstLayer;
    Entry McDirErr;
    Entry McTkr1DirErr;
    Entry CTgoodCal;

    bool m_isThin; 
};

#endif
