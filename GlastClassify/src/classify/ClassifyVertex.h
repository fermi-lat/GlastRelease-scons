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
    {
        m_vtxangleIndex = find_index( "VtxAngle");
        m_firstlayer =     add_index( "Tkr1FirstLayer");
        m_direrr =         add_index( "McDirErr");
        m_tkrdirerr =      add_index( "McTkr1DirErr");
    }

    //function to generate good test

    virtual bool isgood()
    {
        return datum(m_vtxangleIndex) > 0 && datum(m_direrr) < datum(m_tkrdirerr);
    }

    virtual bool accept()
    {
        bool thin(datum(m_firstlayer) > 5);
        return thin==m_isThin;
    }

private:
    int m_vtxangleIndex;
    int m_firstlayer;
    int m_direrr;
    int m_tkrdirerr;
    bool m_isThin; 
};

#endif
