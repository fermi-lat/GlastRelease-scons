// $Id$

#ifndef SI_DATA_H
#define SI_DATA_H 1

// (should-be)nested,  private class knows about each hit

#include "geometry/Point.h"
#include "instrument/ModuleId.h"

#include <iostream>
#include <vector>

class GPStime;
class GlastEvent;
class SiController;
class SiDetector;

// std::vector< SiData_Hit >

class SiData_Hit {
 private:
    friend class SiData;
    friend class SiEventData;

 public:
    SiData_Hit ()
      {
      }

  private:
    SiData_Hit (Point p, ModuleId m, unsigned int n, unsigned int t = 0)
      : stripIndex(n), stripType(t), pos(p), module(m)
      {
      }

      unsigned int stripIndex;
      unsigned int stripType; 
      Point pos;
      ModuleId module;
};


#ifdef _MSC_VER
# pragma warning(disable:4786)
bool operator<(const SiData_Hit&, const SiData_Hit&);
bool operator==(const SiData_Hit&, const SiData_Hit&);
#endif

class SiData
{
public:

    enum Axis{X,Y};

  public:
      SiData ();

      SiData (unsigned int n);

      virtual ~SiData ();

      //	access to module id if desired
      virtual const ModuleId& moduleId (enum SiData::Axis a, 
					unsigned int tray, 
					unsigned int n) const;


      void load (const SiDetector& strips);

      //-------------------- access functions for data analysis------------------------
      void clear ();

      //	number of hit strips in the given tray (tray=0..16)
      int nHits (enum SiData::Axis a, int tray) const;

      //	center of the given strip, where  n=0..nHits(a,tray,n)-1
      Point hit (enum SiData::Axis a, unsigned int tray, unsigned int n) const;

      //	strip Id
      unsigned int hitId (enum SiData::Axis a, 
			  unsigned int tray, unsigned int n) const;

      //	Strip Type (1 = noise, 0 = real)
      unsigned int hitType (enum SiData::Axis a, 
			    unsigned int tray, unsigned int n) const;

      //	total number of Si strips hit
      //	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      int totalHits () const;

      //	largest number of hits on an individual readout
      //	controller
      //	(determines readout time)
      //	----- I/O Functions -----------------------
      //	total number of Si strips hit
      int maxControllerHits () const;

      void readData (std::istream& in);

      //	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      //	Data I/O
      void writeData (std::ostream& out);

      //	---------------------------------------------------------
      void printOn (std::ostream& cout = std::cout) const;

  private:

      int m_total_hits;
      int m_controller_max;
      //	storage is a vector (by layer number) of lists of hits
      std::vector< std::vector< class SiData_Hit >* > xhitList;

      //	storage is a vector (by layer number) of lists of hits
      std::vector< std::vector< class SiData_Hit >* > yhitList;
};

#endif
