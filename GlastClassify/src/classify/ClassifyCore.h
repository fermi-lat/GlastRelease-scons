/**@file ClassifyCore.h
@brief 

*/

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
    {
        m_firstlayer    = add_index( "Tkr1FirstLayer");
        m_direrr        = add_index( "McDirErr");
        m_tkrdirerr     = add_index( "McTkr1DirErr");
        m_energy        = add_index( "EvtEnergyCorr");
        m_CTgoodCal     = add_index( "CTgoodCal");
        m_CTvertex      = add_index( "CTvertex");
    }

    //function to generate good test

    virtual bool isgood()
    {
        // test direrr or tkrdirerr against criterion. adjust the factor to give < 1/2 in "tail"
        return direrror() < 2* scaleFactor(energy(), m_isThin);

    }

    virtual bool accept()
    {
        // test for good energy etc first: goodCal includes this 
        double goodCal =datum(m_CTgoodCal);
        if( datum(m_CTgoodCal)< 0.25) return false;

        // now do this if it is in the correct class
        return ( m_isThin== isThin()) && (useVertex()==m_isVertex);
    }

    /// function that 
    static double scaleFactor(double energy, bool thin);

    double energy()const{return datum(m_energy);}

    /// check to see if vertex is in thin or thick (front or back) sections
    bool isThin() const { return datum(m_firstlayer) > 5;}

    /// true if CT for vertex says to use it.
    bool useVertex() const { return datum(m_CTvertex)>0.2; }


    ///
    double direrror() const{ return useVertex()? datum(m_direrr) : datum(m_tkrdirerr);} 

private:
    int m_firstlayer;
    int m_direrr;
    int m_tkrdirerr;
    int m_energy;
    bool m_isVertex;
    bool m_isThin;
    int m_CTgoodCal;
    int m_CTvertex;
};
