#include "EbfWriter/Ebf.h"

using namespace EbfWriterTds;

Ebf::Ebf(){
    m_data=NULL;
    m_length=0;
}

Ebf::Ebf(char *newData, unsigned int dataLength){
    m_data=NULL;
    m_length=0;
    set(newData,dataLength);
}

Ebf::~Ebf(){
    if(m_data!=NULL)
        delete[] m_data;
}

char *Ebf::get(unsigned int &dataLength) const{
    dataLength=m_length;
    return m_data;
}

inline std::ostream& Ebf::fillStream( std::ostream &s) const{
    if(m_length>0)
        s.write(m_data,m_length);
    return s;
}

std::ostream& operator << (std::ostream& s, const Ebf& obj){
    return obj.fillStream(s);
}

void Ebf::set(char *newData,unsigned int dataLength){
    if(m_data!=NULL)
        delete[] m_data;
    m_data=NULL;
    m_data=new char[dataLength];
    memcpy(m_data,newData,dataLength);
    m_length=dataLength;
}
