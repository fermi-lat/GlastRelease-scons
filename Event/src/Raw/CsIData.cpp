


// Implementatio of the CsIData class for the TDS

#include "GlastEvent/Raw/CsIData.h"
#include "instrument/Calorimeter.h"


CsIData::CsIData (int numLayers)
{
    for(int i=0; i<numLayers; i++)  {    
        calorList.push_back((new std::vector<Xtal>));
    }
}



CsIData::~CsIData ()
{
    clear();
    for(unsigned i=0; i<calorList.size(); i++) {
    	delete calorList[i];
    }
}

int CsIData::nHits (unsigned int layer) const
{
    return calorList[layer]->size();
}

float CsIData::energy (unsigned int layer, unsigned int n) const
{
    return (*calorList[layer])[n].energy;
}

Point CsIData::xtalPos (unsigned int layer, unsigned int n) const
{
    return  (*calorList[layer])[n].pos;
}


ModuleId CsIData::moduleId (unsigned int layer, unsigned int n) const
{
    return (*calorList[layer])[n].module;
}

unsigned int CsIData::xtalId (unsigned int layer, unsigned int n) const
{
    return (*calorList[layer])[n].id;
}

float CsIData::Lresp (unsigned int layer, unsigned int n) const
{
    return (*calorList[layer])[n].Lresp;
}

float CsIData::Rresp (unsigned int layer, unsigned int n) const
{
    return (*calorList[layer])[n].Rresp;
}

/*
const std::vector<double>& CsIData::Diodes_Energy (unsigned int layer, unsigned int n) const
{
    return (*calorList[layer])[n].Diodes_Energy;
}
*/

void CsIData::load (const CsIDetector& xtal)
{
    // The moduleId was taken out and replaced with a zero for now    
    if (xtal.hit()) {
        calorList[xtal.layer()]->push_back(Xtal(xtal.xtalPos(), xtal.energy(),0 ,
            xtal.xtalId(), xtal.Lresp(), xtal.Rresp(),
				xtal.getDiodesEnergy()));

        // put the "raw" information into the basic list
        // Replaced moduleId with a 0 for now must fix
        XtalId id(0, xtal.layer(), xtal.xtalId());
        int	left =  static_cast<int>(xtal.Lresp()*1e6),
            right = static_cast<int>(xtal.Rresp()*1e6);
        m_xtals[id]=std::pair<int,int>(left,right);
    }
}




void CsIData::readData (std::istream& in)
{
   float Left, Right;
   int X, Y, Z, E;

   int numLayers;
   in>> numLayers;

   for(int i=0; i<numLayers; i++) {
       int numX;
       in>>numX;
       unsigned st;
       ModuleId mod;
       for(int j=0; j<numX; j++) {
			in>>mod>>st>>E>>X>>Y>>Z>>Left>>Right;
           calorList[i]->push_back(Xtal(Point(X/1e3, Y/1e3, Z/1e3), E/1e6, mod, st,
			   Left/1e6, Right/1e6));
       }
   }
}

void CsIData::writeData (std::ostream& out) const
{
   int numLayers = calorList.size();

   out<<numLayers<<'\n';
   if(out.eof() ) return;  // happens to be the first
   for(int i=0; i<numLayers; i++) {
       int numX = nHits(i);
       out<<numX<<'\n';
       for(int j=0; j<numX; j++) {
            out<<moduleId(i, j)<<' '<<xtalId(i, j)<<' '
               <<int(1e6*energy(i,j))       <<' '
               <<int(1e3*xtalPos(i, j).x())<<' '
               <<int(1e3*xtalPos(i, j).y())<<' '
               <<int(1e3*xtalPos(i, j).z())<<' '
			   <<int(1e6*Lresp(i, j))<<' '
			   <<int(1e6*Rresp(i, j))<<'\n';
       }
   }
}

void CsIData::clear ()
{
    m_xtals.clear();
    for(unsigned i=0; i<calorList.size(); i++) {
	calorList[i]->clear();
    }
}


