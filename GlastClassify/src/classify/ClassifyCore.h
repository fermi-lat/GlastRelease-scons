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

    ClassifyCore(const std::string& info_path)
        : GlastClassify (info_path)
        , m_isVertex    (info_path.find("vertex")!=std::string::npos)
        , m_isThin      (info_path.find("thin")!=std::string::npos)

        , Tkr1FirstLayer("Tkr1FirstLayer")
        , McDirErr      ("McDirErr")
        , McTkr1DirErr  ("McTkr1DirErr")
        , EvtEnergyCorr ("EvtEnergyCorr")
        , CTgoodCal     ("CTgoodCal")
        , CTvertex      ("CTvertex")

    {
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
        if( CTgoodCal< 0.25) return false;

        // now do this if it is in the correct class
        return ( m_isThin== isThin()) && (useVertex()==m_isVertex);
    }

    /// function that 
    static double scaleFactor(double energy, bool thin);

    double energy()const{return EvtEnergyCorr;}

    /// check to see if vertex is in thin or thick (front or back) sections
    bool isThin() const { return Tkr1FirstLayer > 5;}

    /// true if CT for vertex says to use it.
    bool useVertex() const { return CTvertex>0.2; }


    ///
    double direrror() const{ return useVertex()? McDirErr : McTkr1DirErr;} 

private:
    bool m_isVertex;
    bool m_isThin;
    Entry Tkr1FirstLayer;
    Entry McDirErr;
    Entry McTkr1DirErr;
    Entry EvtEnergyCorr;
    Entry CTgoodCal;
    Entry CTvertex;
};
