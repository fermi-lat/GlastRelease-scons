
//Old School Version

#include "TdSiData.h"

#include "instrument/SiTracker.h"
#include <algorithm>
using namespace std;

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

void TdSiData::load (const SiDetector& plane, idents::ModuleId moduleId)
{
    static double tray_max=48.10, tray_spacing=3.20;
    
    CoordTransform T 
        = plane.localToGlobal(); // for converting the coords following
    
    // figure out which tray this is from the transform
    Vector a(1,0,0),b(1,0,0); a.transform(T);
    double rot = a*b; 
    Axis axis = (rot>0.5)? X : Y;
    Point p(0,0,0); p.transform(T); 
    double z = p.z();
    int tray = static_cast<int>((tray_max-z)/tray_spacing);
    tray = tray< 0?  0: tray;
    tray = tray>17? 17: tray;
    
    // final loop over hits
    for( SiDetector::const_iterator it3=plane.begin(); 
    it3!=plane.end(); ++it3)                      {
        int wire_num = (*it3).index();
        int wire_typ = (*it3).noise();
        //************** remove noise  *********************
        if( wire_typ !=0 ) continue; 
        // *************************************************
        Point p = Point(SiDetector::localX(wire_num),0,0);
        p.transform(T);
        // Right here we are running off the end of the xhitList
        // Toby do you know why
        if( axis == X) { 
            xhitList[tray]->push_back(
                Strip(p, moduleId, wire_num, wire_typ));
        } else {
            yhitList[tray]->push_back(
                Strip(p, moduleId, wire_num, wire_typ));
        }
        m_total_hits++; 
    }
}

void TdSiData::clear () {
    for(unsigned i=0; i<xhitList.size(); i++) {
        xhitList[i]->clear();
        yhitList[i]->clear();
    }
    m_total_hits=0;
}

int TdSiData::nTrays (enum TdSiData::Axis a) const
{
    return a == X ? xhitList.size() : yhitList.size() ;
    
}

int TdSiData::nHits (enum TdSiData::Axis a, int tray) const
{
    int numHits = xhitList[tray]->size();

    if (a != X) numHits = yhitList[tray]->size();

    return numHits;
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
