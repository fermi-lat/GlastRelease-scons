#ifndef CALGAINCALIB_H
#define CALGAINCALIB_H


#include <fstream>

class CalGainElement
{
    public:
        CalGainElement():m_gain(0){ }
        ~CalGainElement() {}
        void generateCalib(const int range = 0)
        { const float adcMidRange = 2000; m_gain = 100.0*pow(8,range)/adcMidRange;}
            
            
        void readFile(std::ifstream& file) { file >> m_gain;}
        void writeFile(std::ofstream& file) const { file << " " << m_gain << std::endl;} 

        const float getGain() const {return m_gain;}
    private:
        float m_gain;
};

class CalGainCalib
{
    public:

        CalGainCalib(){m_elem = (CalGainElement*)0;}
        ~CalGainCalib() {if(m_elem != 0) delete m_elem;}
        void addElement(CalGainElement* elem)
        {
            if(m_elem != 0) delete m_elem;
        
            m_elem = elem;
        } 

        
        void generateCalib(unsigned int range =0){if(m_elem == 0) m_elem = new CalGainElement();
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

        const float getGain() const { return m_elem->getGain();}

    private:
        CalGainElement* m_elem;


};

#endif // CALGAINCALIB_H 
