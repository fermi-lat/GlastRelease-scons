
#ifndef VolumeIdentifier_h
#define VolumeIdentifier_h

#include <string>

/** 
 * @class VolumeIdentifier
 *
 * @brief This class encapsulates volume identifiers.
 *
 * This class encapsulates volume identifiers defined in the xml files 
 * describing the detector. The volume identifiers are transformed into a 
 * 64 bit integer for efficiency purpose.
 *  
 * @author Toby Burnett
 * @author Marco Frailis
 *    
 * \$Header$
 */

namespace idents{


class VolumeIdentifier {
public:

#ifdef WIN32
    typedef __int64 int64;
#else
    typedef  long long int64;
#endif

    VolumeIdentifier();

    /** 
     * This method initialize the VolumeIdentifier with a 64 bit integer 
     * and a size and is used by the overloaded operator<< to build an 
     * identifier reading from a persistent data store
     */
    void init(int64 , unsigned int); 

    int64 getValue() const{ return m_value;}

    /// prepend, in front, another id
    void prepend( const VolumeIdentifier& id);

    /// append another id
    void append( const VolumeIdentifier& id);

    /// append an int
    void append(unsigned int id);

    /// number of single ids which constitute the volume identifier
    int size() const { return m_size;}

    /// return a name made up with slash delimiters
    std::string name(const char * delimiter="/")const ;
 
    /// access single ids which constitute the volume identifier
    unsigned int operator[](unsigned int);
    
    /// access single ids which constitute the volume identifier
    unsigned int operator[](unsigned int) const;
    
    /// overload the < operator for correct sorting of volume identifiers
    bool operator<(const VolumeIdentifier& id)const
      {
         return (m_value < id.getValue())
	   || (m_value==id.getValue() && (m_size < id.size()));
      }
  
    bool operator==(const VolumeIdentifier& id)const
    {
        return ((m_value == id.getValue()) && (m_size==id.size()));
    }

    /// return true iff VolumeIdentifier fields say "tracker"
    bool isTkr() const {return ( (m_size > (int) fTowerObj) &&
                           ( (*this)[fLATObj] == eLATTowers) &&
                           ( (*this)[fTowerObj] == eTowerTKR)  ); }


    /// return true iff VolumeIdentifier fields say "calorimeter"
    bool isCal() const {return ( (m_size > (int) fTowerObj) &&
                           ( (*this)[fLATObj] == eLATTowers) &&
                           ( (*this)[fTowerObj] == eTowerCAL)  ); }

    /// return true iff VolumeIdentifier fields say "ACD"
    bool isAcd() const {return ((m_size > (int) fLATObj) && 
                          ((*this)[fLATObj] == eLATACD)); }

    /// Max allowed value for a single field
    static unsigned maxFieldValue() { return s_maxFieldValue;}
                                                       
private:

    /// The following values must correspond with those defined 
    /// in the xml geometry files in use when the VolumeIdentifier 
    /// was created, or the is.. routines above will lie.
    static const unsigned fLATObj = 0;
    static const unsigned fTowerObj = 3;
    static const unsigned eLATTowers = 0;
    static const unsigned eLATACD   = 1;
    static const unsigned eTowerCAL = 0;
    static const unsigned eTowerTKR = 1;

    static const unsigned s_bitsPer = 6;

    // Since 64 is not evenly divisible by bitsPer (6)  field
    // occupying most significant bits is located in bits 54-59.  
    // Top 4 bits are unused.
    static const unsigned s_maxSize = 64 / s_bitsPer;
    static const unsigned s_maxShift = (s_maxSize - 1) * s_bitsPer;   /* 54 */
    static const unsigned s_maxFieldValue = (1 << s_bitsPer) - 1;

    /// internal rappresentation of the volume identifier
    int64 m_value; // for sorting
    /// number of single ids which constitute the volume identifier
    int m_size;
};

}
#endif
