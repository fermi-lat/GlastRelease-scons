/** 
* @file SpectrumFactoryLoader.cxx
* @brief Load the external spectrum factory objects
*
*  $Header$
*/

#include "SpectrumFactoryLoader.h"
#include "flux/ISpectrumFactory.h"

// declare the external factories.
ISpectrumFactory & GRBmanagerFactory();
ISpectrumFactory & GRBobsFactory();
ISpectrumFactory & IsotropicFactory();
ISpectrumFactory & PulsarSpectrumFactory();


SpectrumFactoryLoader::SpectrumFactoryLoader()
{
    ISpectrumFactory& f =GRBmanagerFactory();
    // fix these in next release, when GRB gets rebuild with new interface
//    m_names.push_back(f.name());
    f=GRBobsFactory();
//    m_names.push_back(f.name());
}
