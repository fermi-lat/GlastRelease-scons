
#include "configData/db/LatcDBImplFile.h"
#include <vector>
#include <iostream>


LatcDBImplFile::LatcDBImplFile():m_filename(""){ }

const std::string LatcDBImplFile::getFilename(const char* type, unsigned key){
  std::cout<<m_filename<<std::endl;
  return m_filename;
}

void LatcDBImplFile::setFilename(const std::string filename){
  m_filename=filename;
}

