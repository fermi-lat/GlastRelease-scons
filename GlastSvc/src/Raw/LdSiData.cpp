// LdSiData.cpp: implementation of the LdSiData class.
//
//////////////////////////////////////////////////////////////////////

#include "LdSiData.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

void LdSiData::load (const SiDetector& plane, idents::ModuleId moduleId)
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
