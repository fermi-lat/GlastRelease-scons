// TdVetoData.cpp: implementation of the TdVetoData class.
//
//////////////////////////////////////////////////////////////////////

#include "TdVetoData.h"
#include "instrument/Scintillator.h"


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


void TdVetoData::clear ()
{
    tileList.clear();
}

