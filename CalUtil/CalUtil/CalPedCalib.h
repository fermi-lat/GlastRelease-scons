#ifndef CALPEDCALIB_H
#define CALPEDCALIB_H


#include <fstream>

class CalPedElement
{
    public:
        CalPedElement():m_pedAvr(0),m_pedSig(0){ }
        ~CalPedElement() {}
        void generateCalib(unsigned int range = 0){m_pedAvr=100.0; m_pedSig=5.0;}
        void readFile(std::ifstream& file) { file >> m_pedAvr >> m_pedSig;}
        void writeFile(std::ofstream& file) const { file << " " << m_pedAvr << " " 
                                             << m_pedSig << std::endl;} 

        const float getAvr() const {return m_pedAvr;}
        const float getSig() const {return m_pedSig;}
    private:
        float m_pedAvr;
        float m_pedSig;
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

        
        void generateCalib(const unsigned int range){if(m_elem == 0) m_elem = new CalPedElement();
                                m_elem->generateCalib(range);}
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
