// TdVetoData.cpp: implementation of the TdVetoData class.
//
//////////////////////////////////////////////////////////////////////

#include "TdVetoData.h"
#include "instrument/Scintillator.h"

//#include "instrument/ScintillatorId.h"

using namespace std;
#include <iomanip>

TdVetoData::TdVetoData ()
{
}


TdVetoData::~TdVetoData ()
{
    clear();
}

void TdVetoData::printOn (ostream& cout) const
{
    cout << "\nTdVetoData:\n";
    for( const_iterator it = begin(); it != end(); ++it){
        const IVetoData::Tile& v = *it;
        cout << '\t' << v.type() << '\t'  << setprecision(3)
            << v.position() << '\t' << v.energy() <<'\n';
    }
}

void TdVetoData::load (const Scintillator& tile)
{
    float energy = tile.energy();
    if( energy >0.0 ) {
        // found a hit one
        CoordTransform T = tile.localToGlobal(); // for converting the coords following
        Point p(0,0,0);
        p.transform(T);
        tileList.push_back(Tile(p, tile.type(), energy) );
    }
    
}


void TdVetoData::readData (istream& in)
{
    int X, Y, Z, E;
    
    int numV;
    Id type;
    in>> numV ;
    
    for(int i=0; i<numV; i++) {
        //      ModuleId mod;
        in >> type >> E >> X >> Y >> Z;
        tileList.push_back( Tile(Point(X/1e3, Y/1e3, Z/1e3), type, E/1e6));
    }
}

void TdVetoData::writeData (ostream& out)
{
    int numV = tileList.size();
    
    out<<numV<<'\n';
    
    for( const_iterator it = begin(); it != end(); ++it) {
        const Tile& v = *it;
        out << v.type()          << ' '
            << static_cast<int>(1e6*v.energy())     << ' '
            << static_cast<int>(1e3*v.position().x())<< ' '
            << static_cast<int>(1e3*v.position().y())<< ' '
            << static_cast<int>(1e3*v.position().z())<< '\n';
    }
}

void TdVetoData::clear ()
{
    tileList.clear();
}

