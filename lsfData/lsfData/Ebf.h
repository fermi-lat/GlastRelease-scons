#ifndef lsfData_Ebf_H
#define lsfData_Ebf_H

#include <iostream>

/**
 * @class Ebf
 *
 * @brief TDS for storing an event in a format similar to ebf
 *
 * The data is stored as one continuos string of bytes
 * No attempt is made to verify that the data stored is correctly
 * formated ebf.
 * $Header$
 */


namespace lsfData {
    class Ebf {
    public:
        Ebf();
        Ebf(char *newData,unsigned int dataLength);
        virtual ~Ebf();

        ///Retrieve pointer to the ebf data.
        char *get(unsigned int &dataLength) const;

        ///Store the provided ebf pointer in and delete any previous ones
        void set(char *newData, unsigned int dataLength);

        unsigned int getSequence() const { return m_gemSeq; };
        void setSequence(unsigned int seq) { m_gemSeq = seq;  };

    private:
        ///Pointer to the ebf data
        char *m_data;
        ///Number of bytes that are stored in data pointer
        unsigned int m_length;
        ///Save the GEM sequence number
        unsigned int m_gemSeq;
    };

    //inline stuff for client
    inline Ebf::Ebf(){ m_data=0; m_length=0;}

    inline  char *Ebf::get(unsigned int &dataLength) const{
      dataLength=m_length;
      return m_data;
    }

    inline Ebf::Ebf(char *newData, unsigned int dataLength){
      m_data=NULL;
      m_length=0;
      set(newData,dataLength);
    }

    inline Ebf::~Ebf(){
      if(m_data!=NULL)
        delete[] m_data;
    }


    inline void Ebf::set(char *newData,unsigned int dataLength){
      if(m_data!=NULL)
        delete[] m_data;
      m_data=NULL;
      m_data=new char[dataLength];
      memcpy(m_data,newData,dataLength);
      m_length=dataLength;
    }
}// namespace
#endif

