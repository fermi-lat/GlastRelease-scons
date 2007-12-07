
#include "configData/db/LatcDBImplFile.h"
#include <vector>
#include <iostream>


LatcDBImplFile::LatcDBImplFile():m_filename(""){ }

const std::string LatcDBImplFile::getFilename(const char* type, unsigned key){
  //  std::cout<<m_filename<<std::endl;
  return m_filename;
}

void LatcDBImplFile::setFilename(const std::string filename){
  m_filename=filename;
}
const std::vector<std::string> LatcDBImplFile::getFilenameList(unsigned key){
  std::vector<std::string> rv;
  if (m_filename!="")rv.push_back(m_filename);
  return rv;
}
