// $Id$

#ifndef VETO_DATA_H
#define VETO_DATA_H 1

#include "geometry/Point.h"
#include "instrument/ScintillatorId.h"

#include <iostream>
#include <iomanip>
#include <vector>

class GPStime;
class GlastEvent;
class VetoController;
class VetoTileData;
class VetoTile;
class Scintillator;
//------------------------------------------------------------------------------
// (not) Nested class holds a single tile

class VetoTileData
{
  public:
      VetoTileData()
      {
      }

      VetoTileData (Point p, ScintillatorId t, float e)
	: m_type(t), m_energy(e), m_pos(p)
      {
      }


      void printOn (std::ostream& os = std::cout) const;

      Point position () const
      {
	return m_pos;
      }

      float energy () const
      {
	return m_energy;
      }

      const ScintillatorId& type () const
      {
	return m_type;
      }


  protected:
  private:
    // Data Members for Class Attributes

      //	type within module
      ScintillatorId m_type;

      float m_energy;

    // Data Members for Associations

      Point m_pos;


  private:
};

class VetoData
{
  public:

    typedef std::vector < VetoTileData >::const_iterator const_iterator;

    class VetoList : public std::vector< VetoTileData >
    {
      public:
      protected:
      private:
      private:
    };


  public:
      VetoData ();

      ~VetoData ();


      void printOn (std::ostream& cout) const;


      //	GlastData interface
      void load (const Scintillator& tile);


      void readData (std::istream& in);

      //	Data I/O
      void writeData (std::ostream& out);

      void clear ();

      VetoData::const_iterator begin () const
      {
	    return tileList.begin();
      }

      VetoData::const_iterator end () const
      {
	    return tileList.end();
      }

  protected:
  private:
    // Data Members for Associations

      //	storage is a list
      VetoList tileList;

  private:
};

#endif
