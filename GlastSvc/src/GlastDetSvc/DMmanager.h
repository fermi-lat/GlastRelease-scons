// $Header$

#ifndef DM_DMmanager_h
#define DM_DMmanager_h

namespace detModel{
    class Manager;
    class Volume;
    class SectionsVisitor;
    class MaterialsVisitor;
    class IDmapBuilder;
}

namespace idents{class VolumeIdentifier;}

class HepTransform3D;

#include <string>

/**
  wrapper class around detModel, to hide all the initialization stuff, provide defaults
  */
class DMmanager  {
public:
    DMmanager();
    ~DMmanager();

    //! initialize with command-line args, supplying defaults
    void init(int argc, char* argv[] ) ;

    //! Note that if filename is "-" or blank, use xmlGeoDbs default
    void init(std::string filename, std::string mode, std::string topvol);

    //! start a visitor from the top volume
    void accept( detModel::SectionsVisitor* v);

    //! start a materials visitor
    void accept( detModel::MaterialsVisitor* v);

    //! retrive a numeric constant
    bool getNumericConstByName(std::string, double*);

    //! summarize setup.
    void printSetup(std::ostream& out);

    /// retrive the 3D transformation of a volume given a valid ID
    bool getTransform3DByID(idents::VolumeIdentifier, HepTransform3D*);

    //! access to name of selected top volume
    std::string topvol()const;

private:
    detModel::Manager* m_dm;
    detModel::Volume * m_vol;
    std::string m_mode;
    /// The idMapBuilder
    detModel::IDmapBuilder* m_idMap;
};

#endif
