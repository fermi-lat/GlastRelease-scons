#ifndef CALPEDCALIB_H
#define CALPEDCALIB_H


#include <fstream>

class CalPedElement
{
    public:
        CalPedElement():pedAvr(0),pedSig(0){ }
        ~CalPedElement() {}
        void generateCalib(){pedAvr=100.0; pedSig=5.0;}
        void readFile(std::ifstream& file) { file >> pedAvr >> pedSig;}
        void writeFile(std::ofstream& file) const { file << " " << pedAvr << " " 
                                             << pedSig << std::endl;} 

        const float getAvr() const {return pedAvr;}
        const float getSig() const {return pedSig;}
    private:
        float pedAvr;
        float pedSig;
};

class CalPedCalib
{
    public:

        CalPedCalib(){m_elem = (CalPedElement*)0;}
        ~CalPedCalib() {if(m_elem != 0) delete m_elem;}
        void addElement(CalPedElement* elem)
        {
            if(m_elem != 0) delete m_elem;
        
            m_elem = elem;
        } 

        
        void generateCalib(){if(m_elem == 0) m_elem = new CalPedElement();
                                m_elem->generateCalib();}
        void writeFile(std::ofstream& file, unsigned int tower,
                                  unsigned int layer,
                                  unsigned int col,
                                  unsigned int range,
                                  idents::CalXtalId::XtalFace face) const
        {
            file << " " << tower
                << " " << layer
                << " " << col
                << " " << range 
                << " " << face;
            m_elem->writeFile(file);
        }

        const float getAvr() const { return m_elem->getAvr();}
        const float getSig() const { return m_elem->getSig();}
    private:
        CalPedElement* m_elem;


};

#endif // CALPEDCALIB_H 
