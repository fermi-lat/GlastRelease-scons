// File and Version Information:
// $Header$
//
// Description: This class implements the abstract interface IMedia of
// GlastDetSvc. The main purpouse is to build the material table for a Geant4
// simulation by reading material specifications from the detModel
// representation (obtained by reading xml files). The use of an abstract
// interface permits to avoid the explicit exposition of detModel interfaces to
// clients. For more detailed information, please see the GlastDetSvc
// documentation
//      
// Author(s): R.Giannitrapani

#include "G4Media.h"

#include "G4Material.hh"
#include "globals.hh"
#include <numeric> // for accumulate


void G4Media::addMaterial(std::string name, 
                          MediaType type,
                          const DoubleVector& params, std::string symbol){
  // Purpose and Method:  This routine add a simple material to the table
  // Inputs: the material name, type (element or material), vector of
  //         parameters (atomic number, density, etc etc) and the symbol 
  //         (in case of elements)
  
  if (type == IMedia::Element)
    {
      new G4Element(name, symbol, params[0], params[1]*g/mole);    
    }
  else 
    {
      // This is temporary waiting for a new material category 
      // with pressure and temp
      if(name==std::string("Vacuum"))
        {
          new G4Material(name, params[0], params[1]*g/mole, 
                         params[2]*g/cm3, kStateGas,
                         2.73*kelvin,3.e-18*pascal);
        }
      else
        new G4Material(name, params[0], params[1]*g/mole, params[2]*g/cm3);
 
    }
}


void G4Media::addComposite(std::string name, 
                           CompositeType type,
                           double density,
                           const StringVector& components, 
                           const DoubleVector& qty){
  // Purpose and Method:  This routine add a composite to the table
  // Inputs: the material name, type (by mass fractions or atomic composition),
  //         density, vector of components and their quantities

  G4Element* ptElement = 0;
  G4Material* ptMaterial = 0;
  
  G4Material* mat = new G4Material(name, density*g/cm3, components.size());
  
  for(unsigned int i=0;i<components.size();i++)
    {
      if (!(ptElement = G4Element::GetElement(components[i])))
        {
          ptMaterial = G4Material::GetMaterial(components[i]);
        }
      
      if (ptElement)
        {
          if (type == IMedia::Fract)
            {
              double tot = std::accumulate(qty.begin(),qty.end(), 0.0);
              mat->AddElement(ptElement, (G4double) qty[i]/tot);
            }
          else
            {
              mat->AddElement(ptElement, (G4int) qty[i]);
            }
        }
      else if (ptMaterial)
        {
          if (type == IMedia::Fract)
            {
              double tot = std::accumulate(qty.begin(),qty.end(), 0.0);
              mat->AddMaterial(ptMaterial, (G4double) qty[i]/tot);
            }
          else
            {
              mat->AddMaterial(ptMaterial, (G4int) qty[i]);
            }
        }
    }
}
