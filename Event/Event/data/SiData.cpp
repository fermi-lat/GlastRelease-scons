// $Id$
//
// Implement access to Silicon strip data

#include "GlastEvent/data/SiData.h"

#include "instrument/SiTracker.h"
#include <algorithm>
using namespace std;

SiData::SiData ()
{
    for(int i = 0; i < SiTracker::numTrays(); i++) {
	xhitList.push_back(new vector<SiData_Hit>);
	yhitList.push_back(new vector<SiData_Hit>);
    }
}

SiData::SiData (unsigned int n)
{
    for(unsigned i=0; i<n; i++) {	
	xhitList.push_back(new vector<SiData_Hit>);
	yhitList.push_back(new vector<SiData_Hit>);
    }
}


SiData::~SiData ()
{
    clear();
    for(unsigned i=0; i<xhitList.size(); i++) {
	delete xhitList[i];
	delete yhitList[i];
    }

}


const ModuleId& SiData::moduleId (enum SiData::Axis a, 
				  unsigned int tray, 
				  unsigned int n      ) const
{
  return a==X? (*xhitList[tray])[n].module
             : (*yhitList[tray])[n].module;
}

void SiData::load (const SiDetector& strips)
{
#if 0
    m_total_hits = 0;
    m_controller_max = 0;
    
    // first loop through all controllers
    for(SiModuleMap::const_iterator it = ctl->begin(); 
	it != ctl->end(); ++it)                        {
      const SiDetectorList& planes = (*it).second;
      ModuleIdXY idxy = (*it).first;
      ModuleId moduleId = idxy.module();
      Axis axis = idxy.axis()==ModuleIdXY::X ? X : Y; // convert to our axis
      int tray = 0; //SiTracker::numTrays(); 
      // assume start at bottom, want tray=0 for top
        
      int controller_hits = 0;
        
      // now loop over planes wired to this controller (from bottom up)
      for( SiDetectorList::const_iterator it2=planes.begin(); 
	   it2!=planes.end(); ++it2)                           {
	const SiDetector& plane = **it2;

	CoordTransform T 
	  = plane.localToGlobal(); // for converting the coords following

	controller_hits += plane.size();

	// final loop over hits
	for( SiDetector::const_iterator it3=plane.begin(); 
	     it3!=plane.end(); ++it3)                      {
	  int wire_num = (*it3).index();
	  int wire_typ = (*it3).noise();
	  Point p = Point(SiDetector::localX(wire_num),0,0);
	  p.transform(T);
	  if( axis == X) { 
	    xhitList[tray]->push_back(
			       SiData_Hit(p, moduleId, wire_num, wire_typ));
	  } else {
	    yhitList[tray]->push_back(
			       SiData_Hit(p, moduleId, wire_num, wire_typ));
	  }
	}
	++tray;
      }
      m_total_hits += controller_hits;
      m_controller_max = std::max(m_controller_max, controller_hits);
    }
#endif
}

void SiData::clear () {
  for(unsigned i=0; i<xhitList.size(); i++) {
    xhitList[i]->clear();
    yhitList[i]->clear();
  }
}

int SiData::nHits (enum SiData::Axis a, int tray) const
{
    return a==X?  xhitList[tray]->size()
		:  yhitList[tray]->size() ;

}

Point SiData::hit (enum SiData::Axis a, 
		   unsigned int tray, 
		   unsigned int n       ) const
{
  return a==X? (*xhitList[tray])[n].pos
             : (*yhitList[tray])[n].pos;
}

unsigned int SiData::hitId (enum SiData::Axis a, 
			    unsigned int tray, 
			    unsigned int n       ) const
{
  return a==X? (*xhitList[tray])[n].stripIndex
             : (*yhitList[tray])[n].stripIndex;
}

unsigned int SiData::hitType (enum SiData::Axis a, 
			      unsigned int tray, 
			      unsigned int n       ) const
{
  return a==X? (*xhitList[tray])[n].stripType
             : (*yhitList[tray])[n].stripType;
}

int SiData::totalHits () const
{
    return m_total_hits;
}

int SiData::maxControllerHits () const
{
     return m_controller_max;
}

void SiData::readData (istream& in)
{
   int x, y, Z;

   int numLayers;
   in>> numLayers;
   if( in.eof() )  // this is the first one
       return;

   for(int i=0; i<numLayers; i++) {
       int numX;
       in>>numX;
       ModuleId mod;
       unsigned st, type;
       int j;
       for( j=0; j<numX; j++) {
           in>>mod>>st>>type>>x>>y>>Z;
           xhitList[i]->push_back(
	       SiData_Hit(Point(x/1e3, y/1e3, Z/1e3), mod, st, type));
       }
       int numY;
       in>>numY;
       for(j=0; j<numY; j++) {
           in>>mod>>st>>type>>x>>y>>Z;
           yhitList[i]->push_back(
	       SiData_Hit(Point(x/1e3, y/1e3, Z/1e3), mod, st, type));
       }
   }
}

void SiData::writeData (ostream& out)
{
   int numLayers = xhitList.size();

   out<<numLayers<<'\n';

   for(int i=0; i<numLayers; i++) {
       int numX = nHits(X, i);
       out<<numX<<'\n';
       int j;
       for(j=0; j<numX; j++) {
	 out<<moduleId(X, i, j)
	    <<' '<<hitId(X, i, j) <<' '<<hitType(X, i, j)<<' '
	    <<int(1e3*hit(X, i, j).x())<<' '
	    <<int(1e3*hit(X, i, j).y())<<' '
	    <<int(1e3*hit(X, i, j).z())<<'\n';
       }
       int numY = nHits(Y, i);
       out<<numY<<'\n';
       for(j=0; j<numY; j++) {
	 out<<moduleId(Y, i, j) 
	    <<' '<<hitId(Y, i, j)<<' '<<hitType(Y, i, j)<<' '
	    <<int(1e3*hit(Y, i, j).x())<<' '
	    <<int(1e3*hit(Y, i, j).y())<<' '
	    <<int(1e3*hit(Y, i, j).z())<<'\n';
       }
   }
}

void SiData::printOn (ostream& cout) const
{
  unsigned tray, i;
  cout << "\nSiData:\n";
  for(tray=0; tray < xhitList.size() && tray < yhitList.size(); tray++)
  {
     unsigned nx= nHits(X,tray),
         ny= nHits(Y,tray);
     if( nx==0 && ny==0 ) continue;
     cout <<  tray << "X";
     for( i =0; i<nx; i++)
        cout << '\t'<< hit(X,tray,i) <<"\t " << hitType(X,tray,i) << '\n';
     cout <<  tray << "Y";
     for( i =0; i<ny; i++)
        cout << '\t'<< hit(Y,tray,i) <<"\t " << hitType(Y,tray,i) << '\n';
  }
  return;
}

// Class SiData_Hit




