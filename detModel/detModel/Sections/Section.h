#ifndef SECTION_H
#define SECTION_H
#include <vector>
#include "detModel/Management/SectionsVisitor.h"

namespace detModel{
  class Choice;
  class Volume;

  /**
   * @author D.Favretto & R.Giannitrapani
   */
  class Section {
  public:
    Section(std::string pname,
	    std::string pversion,
	    std::string pdate,
	    std::string pauthor,
	    std::string ptopVolumeRef):name(pname),version(pversion),date(pdate),
      author(pauthor),topVolumeRef(ptopVolumeRef),
      topVolume(0){};
  
  
    ~Section();

    void addVolume(Volume* nextVolume){volumes.push_back(nextVolume);}
    void addChoice(Choice* nextChoice){choices.push_back(nextChoice);}
    void Accept(SectionsVisitor* v);
    void AcceptNotRec(SectionsVisitor* v);

    std::string getName(){return name;}    
    std::string getVersion()const{return version;};
    std::string getDate()const{return date;};
    std::string getAuthor()const{return author;};
    Volume* getTopVolume()const{return topVolume;};
    std::string getTopVolumeRef()const{return topVolumeRef;};

    void setTopVolume(Volume* pvolume){topVolume = pvolume;};
    void setVersion(std::string pver){version = pver;};
    void setDate(std::string pdate){date = pdate;};
    void setAuthor(std::string pauth){author = pauth;};
    void setTopVolumeRef(std::string pref){topVolumeRef = pref;}

    std::vector< Volume * > getVolumes()const{return volumes;}
    std::vector < Choice * > getChoices()const{return choices;}



  private:
    /** @link aggregation
     *  @supplierCardinality 1..n */
    std::vector < Volume * > volumes;
    /** @link aggregation */
    std::vector < Choice * > choices;
    std::string name;
    std::string version;
    std::string date;
    std::string author;
    std::string DTDVersion;
    std::string topVolumeRef;
    Volume* topVolume;
  };
}
#endif //SECTION_H
