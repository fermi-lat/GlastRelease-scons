// $Id$
//
// Implement access to Silicon strip data

#include "GlastEvent/Raw/TdSiData.h"

#include "instrument/SiTracker.h"
#include <algorithm>
using namespace std;

/*
class TdSiData::Strip {
private:
    friend class TdSiData;
    
public:
    Strip () {}
    
private:
    Strip (Point p, Id m, unsigned int n, unsigned int t = 0)
        : stripIndex(n), stripType(t), pos(p), module(m)
    {}
    
    unsigned int stripIndex;
    unsigned int stripType; 
    Point pos;
    Id module;
};
*/

TdSiData::TdSiData ()
{
    for(int i = 0; i < SiTracker::numTrays(); i++) {
        xhitList.push_back(new vector<Strip>);
        yhitList.push_back(new vector<Strip>);
    }
    m_total_hits=0;
}

TdSiData::TdSiData (unsigned int n)
{
    for(unsigned i=0; i<n; i++) {	
        xhitList.push_back(new vector<Strip>);
        yhitList.push_back(new vector<Strip>);
    }
}

/*! This is the copy constructor. It's used in the converter for the class.
    This is because we need to make a copy of the object which is a member
    of the IRFConverter. Otherwise we lose the object when the destructor
    of IRFConverter is called.
*/
TdSiData::TdSiData ( TdSiData* copy, unsigned int n)
{

}

TdSiData::~TdSiData ()
{
    clear();
    for(unsigned i=0; i<xhitList.size(); i++) {
        delete xhitList[i];
        delete yhitList[i];
    }
    
}


const SiData::Id& TdSiData::moduleId (enum TdSiData::Axis a, 
                                  unsigned int tray, 
                                  unsigned int n      ) const
{
    return a==X? (*xhitList[tray])[n].module
        : (*yhitList[tray])[n].module;
}


void TdSiData::clear () {
    for(unsigned i=0; i<xhitList.size(); i++) {
        xhitList[i]->clear();
        yhitList[i]->clear();
    }
    m_total_hits=0;
}

int TdSiData::nHits (enum TdSiData::Axis a, int tray) const
{
    return a==X?  xhitList[tray]->size()
        :  yhitList[tray]->size() ;
    
}

Point TdSiData::hit (enum TdSiData::Axis a, 
                   unsigned int tray, 
                   unsigned int n       ) const
{
    return a==X? (*xhitList[tray])[n].pos
        : (*yhitList[tray])[n].pos;
}

unsigned int TdSiData::hitId (enum TdSiData::Axis a, 
                            unsigned int tray, 
                            unsigned int n       ) const
{
    return a==X? (*xhitList[tray])[n].stripIndex
        : (*yhitList[tray])[n].stripIndex;
}
#if 0
unsigned int TdSiData::hitType (enum TdSiData::Axis a, 
                              unsigned int tray, 
                              unsigned int n       ) const
{
    return a==X? (*xhitList[tray])[n].stripType
        : (*yhitList[tray])[n].stripType;
}
#endif
int TdSiData::totalHits () const
{
    return m_total_hits;
}

void TdSiData::readData (istream& in)
{
    int x, y, Z;
    
    int numLayers;
    in>> numLayers;
    if( in.eof() )  // this is the first one
        return;
    
    for(int i=0; i<numLayers; i++) {
        int numX;
        in>>numX;
        unsigned int mod;
        unsigned st, type;
        int j;
        for( j=0; j<numX; j++) {
            in>>mod>>st>>type>>x>>y>>Z;
            xhitList[i]->push_back(
                Strip(Point(x/1e3, y/1e3, Z/1e3), idents::ModuleId(mod), st, type));
        }
        int numY;
        in>>numY;
        for(j=0; j<numY; j++) {
            in>>mod>>st>>type>>x>>y>>Z;
            yhitList[i]->push_back(
                Strip(Point(x/1e3, y/1e3, Z/1e3), idents::ModuleId(mod), st, type));
        }
    }
}

void TdSiData::writeData (ostream& out)
{
    int numLayers = xhitList.size();
    
    out<<numLayers<<'\n';
    
    for(int i=0; i<numLayers; i++) {
        int numX = nHits(X, i);
        out<<numX<<'\n';
        int j;
        for(j=0; j<numX; j++) {
            out<<moduleId(X, i, j)
                <<' '<<hitId(X, i, j) <<' '
//                <<hitType(X, i, j)<<' '
                <<int(1e3*hit(X, i, j).x())<<' '
                <<int(1e3*hit(X, i, j).y())<<' '
                <<int(1e3*hit(X, i, j).z())<<'\n';
        }
        int numY = nHits(Y, i);
        out<<numY<<'\n';
        for(j=0; j<numY; j++) {
            out<<moduleId(Y, i, j) 
                <<' '<<hitId(Y, i, j)<<' '
//                <<hitType(Y, i, j)<<' '
                <<int(1e3*hit(Y, i, j).x())<<' '
                <<int(1e3*hit(Y, i, j).y())<<' '
                <<int(1e3*hit(Y, i, j).z())<<'\n';
        }
    }
}

void TdSiData::printOn (ostream& cout) const
{
    unsigned tray, i;
    cout << "\nTdSiData:\n";
    for(tray=0; tray < xhitList.size() && tray < yhitList.size(); tray++)
    {
        unsigned nx= nHits(X,tray),
            ny= nHits(Y,tray);
        if( nx==0 && ny==0 ) continue;
        cout <<  tray << "X";
        for( i =0; i<nx; i++)
            cout << '\t'<< hit(X,tray,i) 
//            <<"\t " << hitType(X,tray,i) 
            << '\n';
        cout <<  tray << "Y";
        for( i =0; i<ny; i++)
            cout << '\t'<< hit(Y,tray,i)  
//            <<"\t " << hitType(Y,tray,i) 
            << '\n';
    }
    return;
}

// Class TdSiData_Hit




