/**@file ClassifyCore.h
@brief 

*/

#pragma once
#include "classifier/Classifier.h"
#include "GlastClassify.h"
#include "ClassifyVertex.h"

#include <cmath>

class ClassifyCore : public GlastClassify
{
public:
	typedef enum {TRACK, VERTEX} etype;

	ClassifyCore(const std::string& info_path, etype isVertex, ClassifyVertex::layer isThin)
        : GlastClassify(info_path)
        , m_isVertex(isVertex)
        , m_isThin(isThin)
    {}

    //functions to check or declare variables

    virtual void define(std::vector<std::string>& all_names)
    {
        m_vtxangleIndex = subdefine(all_names, "VtxAngle");
        m_firstlayer = subdefine(all_names, "Tkr1FirstLayer");
        m_direrr = subdefine(all_names, "McDirErr");
        m_tkrdirerr = subdefine(all_names, "McTkr1DirErr");
        m_energy = subdefine(all_names, "EvtEnergyCorr");
    }

    //function to generate good test

    virtual bool isgood()
    {
        // test direrr or tkrdirerr against criterion. adjust the factor to give < 1/2 in "tail"
        return 3*direrror()<scaleFactor(energy(), m_isThin);

    }

    virtual bool accept()
    {
        // test for good energy etc first      
        bool thin(datum(m_firstlayer) > 5);
		bool usevertex(datum(m_vtxangleIndex)>0.);
        return thin==m_isThin && usevertex==m_isVertex;
    }

    static double scaleFactor(double energy, bool thin);

    double energy()const{return datum(m_energy);}
    double direrror() const{ return m_isVertex? datum(m_direrr) : datum(m_tkrdirerr);} 

private:
    int m_vtxangleIndex;
    int m_firstlayer;
    int m_direrr;
    int m_tkrdirerr;
    int m_energy;
    bool m_isVertex;
    bool m_isThin;
};
