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
	typedef enum {THICK, THIN} layer;

    ClassifyVertex(const std::string& info_path, layer isThin)
        : GlastClassify(info_path)
        , m_isThin(isThin)
        , VtxAngle      ( "VtxAngle")
        , Tkr1FirstLayer( "Tkr1FirstLayer")
        , McDirErr      ( "McDirErr")
        , McTkr1DirErr  ( "McTkr1DirErr")
    {
    }

    //function to generate good test

    virtual bool isgood()
    {
        return VtxAngle > 0 && McDirErr < McTkr1DirErr;
    }

    virtual bool accept()
    {
        bool thin(Tkr1FirstLayer > 5);
        return thin==m_isThin;
    }

private:
    Entry VtxAngle;      
    Entry Tkr1FirstLayer;
    Entry McDirErr;
    Entry McTkr1DirErr;

    bool m_isThin; 
};

#endif
