// $Header$
// 
//  Original author: R.Giannitrapani
//

#include <sstream>
#include <map>
#include <vector>

#include "XmlStreamer.h"

#include "zlib.h"
#include "HepRepSvc/IRegistry.h"
#include "HepRepSvc/IFiller.h"
#include "XmlBuilder.h"

// This method implement the abstract one that save the HepRep to a
// file; in that case to a compressed XML file
void XmlStreamer::saveHepRep(std::string nameFile)
{
  // Note that the extension is set here
  gzFile file = gzopen((nameFile + ".xml.gz").c_str(), "wb5");  

  // Set the builder for the registry @todo I need to change in order
  // to build one builder at the start and use always it
  XMLBuilder builder(file);
  m_registry->useBuilder(&builder);
  
  typedef std::vector<IFiller*> fillerCol;
  std::vector<std::string>::const_iterator i; 

  gzprintf(file,"<heprep>\n");
  
  const std::vector<std::string>& types = m_registry->getTypeTrees();
  
  for(i=types.begin();i!=types.end();i++)
    {
      fillerCol::const_iterator j;
      gzprintf(file,
               "<typetree name=\"%s\" version=\"1.0\">\n",  
               (*i).c_str());
      
      fillerCol temp = m_registry->getFillersByType((*i));
      
      for(j=temp.begin();j!=temp.end();j++)
        {
          (*j)->buildTypes();
          builder.endTypes();
        }

      gzprintf(file,"</typetree>\n");
      builder.reset();
    }

  std::map<std::string, std::string>::const_iterator it;

  std::vector<std::string> typesList;

  for(it=m_registry->getInstanceTrees().begin();
      it!=m_registry->getInstanceTrees().end();
      it++)
    {
      IRegistry::fillerCol::const_iterator jt;
      IRegistry::fillerCol temp = m_registry->getFillersByType(it->second);
      gzprintf(file,
              "<instancetree name=\"%s\" version=\"1.0\" reqtypetree=\"\" ",  
              it->second.c_str()); 
      gzprintf(file,"typeName=\"%s\" typeVersion=\"1.0\">\n" ,it->first.c_str());

      for(jt=temp.begin();jt!=temp.end();jt++)
        {
          (*jt)->fillInstances(typesList);          
          builder.endInstances();
        }
      gzprintf(file,"</instancetree>\n");
      builder.reset();
    }      

  gzprintf(file,"</heprep>\n");
  
  gzclose(file);

}


