// $Header$
//
// Implement access to Veto Scintillator data

#include "GlastEvent/data/VetoData.h"

using namespace std;

VetoData::VetoData ()
{
}


VetoData::~VetoData ()
{
   clear();
}

void VetoData::printOn (ostream& cout) const
{
    cout << "\nVetoData:\n";
    for( const_iterator it = begin(); it != end(); ++it){
	const VetoTileData& v = *it;
	cout << '\t' << v.type() << '\t'  << setprecision(3)
	    << v.position() << '\t' << v.energy() <<'\n';
    }
}

void VetoData::load (const Scintillator& tile)
{
#if 0
    for( VetoTileList::const_iterator it2=ctl->begin(); it2!=ctl->end(); ++it2){
        const Scintillator& tile = **it2;
        float energy = tile.energy();
        if( energy >0.0 ) {
            // found a hit one
            CoordTransform T = tile.localToGlobal(); // for converting the coords following
            Point p(0,0,0);
            p.transform(T);
            tileList.push_back(VetoTileData(p, tile.type(), energy) );
        }
    }
#endif

}




void VetoData::readData (istream& in)
{
    int X, Y, Z, E;

    int numV;
    ScintillatorId type;
    in>> numV ;

    for(int i=0; i<numV; i++) {
      //	ModuleId mod;
        in >> type >> E >> X >> Y >> Z;
        tileList.push_back( VetoTileData(Point(X/1e3, Y/1e3, Z/1e3), type, E/1e6));
    }
}

void VetoData::writeData (ostream& out)
{
   int numV = tileList.size();

    out<<numV<<'\n';

    for( const_iterator it = begin(); it != end(); ++it) {
	const VetoTileData& v = *it;
	out << v.type()          << ' '
	    << static_cast<int>(1e6*v.energy())     << ' '
	    << static_cast<int>(1e3*v.position().x())<< ' '
            << static_cast<int>(1e3*v.position().y())<< ' '
            << static_cast<int>(1e3*v.position().z())<< '\n';
    }
}

void VetoData::clear ()
{
    tileList.clear();
}


void VetoTileData::printOn (ostream& os) const
{
}

