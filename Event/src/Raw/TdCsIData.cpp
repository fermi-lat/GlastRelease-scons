

// Implementatio of the TdCsIData class for the TDS

#include "GlastEvent/Raw/TdCsIData.h"
//#include "instrument/Calorimeter.h"

//! Default constructor
TdCsIData::TdCsIData (int numLayers)
{
    for(int i=0; i<numLayers; i++)  {    
        calorList.push_back((new std::vector<Xtal>));
    }
}

/*! Method to copy  all the relevant information from another TdCsIData
    object.
*/
void TdCsIData::copyUp( TdCsIData* copy,int numLayers)
{
    for(int i = 0; i < numLayers; i++)
    {
        for(int j = 0;j  < copy->nHits(i); j++)
        {
            // Taking out the diodes for now...fix later.
            calorList[i]->push_back(Xtal(copy->xtalPos(i,j), copy->energy(i,j),0/*moduleId*/,
            copy->xtalId(i,j), copy->Lresp(i,j), copy->Rresp(i,j)/*copy->Diodes_Energy(i,j)*/));
        }
    }
}

//! Destructor
TdCsIData::~TdCsIData ()
{
    clear();
    for(unsigned i=0; i<calorList.size(); i++) {
    	delete calorList[i];
    }
}

//! retern the number of Hits of a specific layer
int TdCsIData::nHits (unsigned int layer) const
{
    return calorList[layer]->size();
}

//! renurn energy on specific xtal at layer and index
float TdCsIData::energy (unsigned int layer, unsigned int n) const
{
    if(n < calorList[layer]->size())
    {
        return (*calorList[layer])[n].energy;
    } else
        return -2;
}

//! return position of specific xtal at layer and index
Point TdCsIData::xtalPos (unsigned int layer, unsigned int n) const
{
    return  (*calorList[layer])[n].pos;
}

// Will implement later
/*ModuleId TdCsIData::moduleId (unsigned int layer, unsigned int n) const
{
    return (*calorList[layer])[n].module;
}*/

//! get the id 
unsigned int TdCsIData::xtalId (unsigned int layer, unsigned int n) const
{
    return (*calorList[layer])[n].id;
}

//! get the  Lresp of specific xtal at layer and index
float TdCsIData::Lresp (unsigned int layer, unsigned int n) const
{
    return (*calorList[layer])[n].Lresp;
}

//! get the Rresp of specific xtal at layer and index
float TdCsIData::Rresp (unsigned int layer, unsigned int n) const
{
    return (*calorList[layer])[n].Rresp;
}


const std::vector<double>& TdCsIData::Diodes_Energy (unsigned int layer, unsigned int n) const
{
    return (*calorList[layer])[n].Diodes_Energy;
}

ModuleId TdCsIData::moduleId (unsigned int layer, unsigned int n) const
{
    return (*calorList[layer])[n].module;
}


/*! Takes in a TdCsIData object makes an xtal object and pushes in back onto
    the  calorList of xtal objects. Currently called from IRFConverter
*/
void TdCsIData::load (const CsIDetector& xtal)
{
    
    // The moduleId was taken out and replaced with a zero for now    
    if (xtal.hit()) {

        Xtal* george = new Xtal(xtal.xtalPos(), xtal.energy(),0/*moduleId*/,
            xtal.xtalId(), xtal.Lresp(), xtal.Rresp(),
				xtal.getDiodesEnergy());

        calorList[xtal.layer()]->push_back(Xtal(xtal.xtalPos(), xtal.energy(),0/*moduleId*/,
            xtal.xtalId(), xtal.Lresp(), xtal.Rresp(),
				xtal.getDiodesEnergy()));
        int temp = xtal.layer();
        // put the "raw" information into the basic list
        // Replaced moduleId with a 0 for now must fix
    //       XtalId id(0, xtal.layer(), xtal.xtalId());
        int	left =  static_cast<int>(xtal.Lresp()*1e6),
            right = static_cast<int>(xtal.Rresp()*1e6);
    //        m_xtals[id]=std::pair<int,int>(left,right);
    }


}





void TdCsIData::writeData (std::ostream& out) const
{
   int numLayers = calorList.size();

   out<<numLayers<<'\n';
   if(out.eof() ) return;  // happens to be the first
   for(int i=0; i<numLayers; i++) {
       int numX = nHits(i);
       out<<numX<<'\n';
       for(int j=0; j<numX; j++) {
            out<</*moduleId(i, j)<<*/' '<<xtalId(i, j)<<' '
               <<int(1e6*energy(i,j))       <<' '
               <<int(1e3*xtalPos(i, j).x())<<' '
               <<int(1e3*xtalPos(i, j).y())<<' '
               <<int(1e3*xtalPos(i, j).z())<<' '
			   <<int(1e6*Lresp(i, j))<<' '
			   <<int(1e6*Rresp(i, j))<<'\n';
       }
   }
}

void TdCsIData::clear ()
{
//    m_xtals.clear();
    for(unsigned i=0; i<calorList.size(); i++) {
	calorList[i]->clear();
    }
}



