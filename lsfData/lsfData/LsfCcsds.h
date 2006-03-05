#ifndef lsfData_LsfCcsds_H
#define lsfData_LsfCcsds_H

namespace lsfData {

    /** @class LsfCcsds
      * @brief Local storage of error data
      * $Header$
    */
    class LsfCcsds {
    public:

        LsfCcsds() { clear(); };
        LsfCcsds(const LsfCcsds& ccsds) { 
            m_scid = ccsds.m_scid; 
            m_apid = ccsds.m_apid; 
            m_utc = ccsds.m_utc; };

        ~LsfCcsds() { clear(); };

        void clear() { 
            m_scid = 0;
            m_apid = 0;
            m_utc = 0.0;
        };

       void print() const {
        }
      
       void initialize(int scid, int apid, double utc) {
           m_scid = scid;
           m_apid = apid;
           m_utc = utc;
       }

       int getScid() const { return m_scid; };
 
       int getApid() const { return m_apid; };
   
       double getUtc() const { return m_utc; };


    private:

        int m_scid;
        int m_apid;
        double m_utc;
    };
} // end namespace
#endif
