#ifndef ACDPOCASORTER_H
#define ACDPOCASORTER_H

#include "../AcdRecon/AcdReconStruct.h"

class AcdPocaSorter {

public:
  enum PocaType { PlanePoca = 0, RayPoca = 1 };
  enum TrackDirection { Upward = 0, Downward = 1 };

  class AcdPocaHolder {    
  public:
    AcdPocaHolder(PocaType type, AcdRecon::PocaData* poca)
      :m_type(type),m_poca(poca){;}
    AcdPocaHolder(const AcdPocaHolder& other)
      :m_type(other.m_type),m_poca(other.m_poca){;}      
      
    ~AcdPocaHolder(){;}
    // access
    inline PocaType type() const { return m_type; }
    inline AcdRecon::PocaData* poca() const { return m_poca; } 
    double arclength() const;
    
    bool operator==(const AcdPocaHolder& other) const;
    bool operator<(const AcdPocaHolder& other) const;
    bool operator>(const AcdPocaHolder& other) const;
    bool operator<=(const AcdPocaHolder& other) const;
    bool operator>=(const AcdPocaHolder& other) const;   

  private:  
    PocaType m_type;
    AcdRecon::PocaData* m_poca;
  };
  
  AcdPocaSorter(TrackDirection dir, const AcdRecon::PocaDataPtrMap& theMap); 
  ~AcdPocaSorter(){;}
  unsigned getPocasToArclength(const double& stop, std::vector<AcdPocaSorter::AcdPocaHolder>& pocas);
  void resetArcCache();

protected:
  void sort();
  void addPoca(const AcdRecon::PocaData* pocaData);
  
private:

  std::vector<AcdPocaSorter::AcdPocaHolder>::iterator m_cache;
  std::vector<AcdPocaSorter::AcdPocaHolder>::reverse_iterator m_rcache;
  std::vector<AcdPocaSorter::AcdPocaHolder> m_pocas;
  TrackDirection m_dir;
  
};




#endif
