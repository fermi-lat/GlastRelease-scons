
#ifndef reconRootData_CalClassParams_H
#define reconRootData_CalClassParams_H

#include <TObject.h>
#include <map>

/** 
* @class CalClassParams
*
* @brief Root clone for the TDS CalClassParams
* @author Luca Baldini, Johan Bregeon
*
*/

class CalClassParams
{
public:

  /// Default constructor
  CalClassParams();

  /// Constructor from all members.
  CalClassParams(std::string producerName, std::map <std::string, double> probMap);
  
  /// Destructor.
  virtual ~CalClassParams() {}
  
  /// Initialize the class members.
  void init(std::string producerName, std::map <std::string, double> probMap);

  inline std::string getProducerName()                 const { return m_producerName; }
  inline std::map <std::string, Double_t> getProbMap() const { return m_probMap; }

  inline void setProbMap(std::map <std::string, double> probMap) { m_probMap = probMap; }
  inline void setProducerName(std::string producerName)      { m_producerName = producerName; }
  inline void setProb(std::string className, double prob)    { m_probMap[className] = prob; }

  Double_t getProb(const std::string &className)       const;
  Double_t getGamProb()                                const { return getProb("gam") ; }
  Bool_t hasClass(const std::string &className)        const;
  
  void Clear( Option_t * option ="" ) ;
  void Fake( Int_t ievent, UInt_t rank, Float_t randNum ) ; // for tests
  Bool_t CompareInRange( const CalClassParams &, const std::string & name = "" ) const ; // for tests
  void Print( Option_t * option ="" ) const ;
  
 private:
  
  std::string m_producerName;
  std::map <std::string, double> m_probMap;
  
  ClassDef(CalClassParams,1)
} ;




#endif

