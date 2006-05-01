/*!
*
* \brief
*     Class for mapping an from electronics space to detector space and back.
*
* \author
*     E. Charles, from  ACDmap
*       R. Claus -- GLAST LAT I&T/Online - (claus@slac.stanford.edu)
*
* $Revision$
* \date   February 24, 2004 -- Created
*
*     $Id$
*
* Copyright:
*                               Copyright 2006
*                                     by
*                        The Board of Trustees of the
*                      Leland Stanford Junior University.
*                             All rights reserved.
*
*/

#ifndef AcdMap_H
#define AcdMap_H

#include <map>

namespace AcdMap {

  /*! \class AcdTilePmt
   *  \brief Helper class for accessing the items in the AcdMap look-up table
   */
  
  class AcdTilePmt {
  public:
    inline AcdTilePmt()
      :_pmt(0xFFFF),_tile(0xFFFF){
    }
    inline AcdTilePmt(unsigned p, unsigned t) 
      :_pmt(p),_tile(t){}
    inline ~AcdTilePmt(){}
    
    inline unsigned    tile() const { return _tile; }
    inline unsigned    pmt()  const { return _pmt; }
    
    inline bool operator==(const AcdTilePmt& other) const {
      return tile() == other.tile() && pmt() == other.pmt();
    }
    inline bool operator < (const AcdTilePmt& other) const {
      if ( tile() < other.tile() ) return true;
      if ( tile() > other.tile() ) return false;
      if ( pmt() < other.pmt() ) return true;
      return false;
    }
    inline bool operator > (const AcdTilePmt& other) const {
      if ( tile() > other.tile() ) return true;
      if ( tile() < other.tile() ) return false;
      if ( pmt() > other.pmt() ) return true;
      return false;
    }

  private:
    unsigned    _pmt;
    unsigned    _tile;
  };

  
  /*! \class AcdGarcGafe
   *  \brief Helper class for accessing the items in the AcdMap look-up table
   */
  
  class AcdGarcGafe {
  public:
    inline AcdGarcGafe()
      :_garc(0xFFFF),_gafe(0xFFFF){ 
    }
    inline AcdGarcGafe(unsigned index)
      :_garc(index/18),_gafe(index%18){    
    }
    inline AcdGarcGafe(unsigned garc, unsigned gafe)
      :_garc(garc),_gafe(gafe){
    }
    inline ~AcdGarcGafe() {
    }
    inline unsigned garc() const { return _garc; }
    inline unsigned gafe() const { return _gafe; }
    inline unsigned index() const { return _gafe + 18*_garc; }
    
    inline bool operator==(const AcdGarcGafe& other) const {
      return index() == other.index();
    }
    inline bool operator<(const AcdGarcGafe& other) const {
      return index() < other.index();
    }
    inline bool operator>(const AcdGarcGafe& other) const {
      return index() > other.index();
    }
  private:
    unsigned _garc;
    unsigned _gafe;
  };


  inline AcdTilePmt* _initTbl();
  inline unsigned* _initGemTbl();
  const std::map<AcdTilePmt,AcdGarcGafe> _initMap();
  const std::map<unsigned,unsigned> _initGemMap();

  const AcdTilePmt* toTilePmt(unsigned  cable, unsigned  channel) {
    if (channel >= 18 || cable >= 12 ) return 0;
    AcdTilePmt* table = _initTbl();
    return &(table[cable*18 + channel]);
  };

  void convertToTilePmt(unsigned  cable, unsigned  channel, unsigned& tile, unsigned& pmt) {
    const AcdTilePmt* tilePmt = toTilePmt(cable,channel);
    if ( tilePmt == 0 ) {
      tile = 0xFFFF; pmt = 0xFFFF;	
    } else {
      tile = tilePmt->tile(); pmt = tilePmt->pmt();
    }
  }

  const AcdGarcGafe* toGarcGafe(unsigned tile, unsigned pmt) {
    const std::map<AcdTilePmt,AcdGarcGafe>& map = _initMap();
    std::map<AcdTilePmt,AcdGarcGafe>::const_iterator itr = map.find( AcdTilePmt(pmt,tile) );
    if ( itr == map.end() ) return 0;
    return &(itr->second);
  };

  void convertToGarcGafe(unsigned tile, unsigned pmt, unsigned& garc, unsigned& gafe) {
    const AcdGarcGafe* gg = toGarcGafe(tile,pmt);
    if ( gg == 0 ) {
      garc = 0xFFFF; gafe = 0xFFFF;
    } else {
      garc = gg->garc(); gafe = gg->gafe();
    }
  }
  

  unsigned           tileFromGemIndex(unsigned gemIndex) {
    if ( gemIndex >= 128 ) return 0xFFFF;
    const unsigned* table = _initGemTbl();
    return table[gemIndex];
  };

  unsigned           gemIndexFromTile(unsigned tile) {
    const std::map<unsigned,unsigned>& map = _initGemMap();
    std::map<unsigned,unsigned>::const_iterator itr = map.find(tile);
    if ( itr == map.end() ) return 0xFFFF;
    return itr->second;
  };

}



#ifndef __CINT__ 

namespace AcdMap {

  /*!
   * \Brief Private` function for initializing the ACD map.
   *
   *  This function really should be a private method in ACDmap but Visual Studio
   *  gives an Internal Compiler Error when compiling LDFdump.cpp.
   *
   */
  inline AcdTilePmt* _initTbl() {
    static AcdTilePmt a[12][18] =
      {{AcdTilePmt(0, 120),  AcdTilePmt(0, 1004), // 0
	AcdTilePmt(0, 502),  AcdTilePmt(0, 1005),
	AcdTilePmt(0, 121),  AcdTilePmt(0, 111),
	AcdTilePmt(0, 101),  AcdTilePmt(0, 122),
	AcdTilePmt(0, 112),  AcdTilePmt(0, 102),
	AcdTilePmt(0, 103),  AcdTilePmt(0, 113),
	AcdTilePmt(0, 123),  AcdTilePmt(0, 503),
	AcdTilePmt(0, 104),  AcdTilePmt(0, 114),
	AcdTilePmt(0, 124),  AcdTilePmt(0, 130)},
       {AcdTilePmt(1, 130),  AcdTilePmt(1, 120), // 1
	AcdTilePmt(1, 110),  AcdTilePmt(1, 100),
	AcdTilePmt(1, 500),  AcdTilePmt(1, 121),
	AcdTilePmt(1, 111),  AcdTilePmt(1, 101),
	AcdTilePmt(1, 102),  AcdTilePmt(1, 112),
	AcdTilePmt(1, 122),  AcdTilePmt(1, 103),
	AcdTilePmt(1, 113),  AcdTilePmt(1, 123),
	AcdTilePmt(1, 1005),  AcdTilePmt(1, 501),
	AcdTilePmt(1, 1004),  AcdTilePmt(1, 124)},
       {AcdTilePmt(0, 110),  AcdTilePmt(0, 100), // 2
	AcdTilePmt(0, 220),  AcdTilePmt(0, 210),
	AcdTilePmt(0, 600),  AcdTilePmt(0, 200),
	AcdTilePmt(0,   0),  AcdTilePmt(0,  10),
	AcdTilePmt(0,  20),  AcdTilePmt(0, 221),
	AcdTilePmt(0, 211),  AcdTilePmt(0, 201),
	AcdTilePmt(0,   1),  AcdTilePmt(0,  11),
	AcdTilePmt(0,  21),  AcdTilePmt(0,  22),
	AcdTilePmt(0, 1006),  AcdTilePmt(0,   2)},
       {AcdTilePmt(1,   2),  AcdTilePmt(1, 601), // 3
	AcdTilePmt(1,  22),  AcdTilePmt(1,  21),
	AcdTilePmt(1,  11),  AcdTilePmt(1,   1),
	AcdTilePmt(1, 201),  AcdTilePmt(1, 211),
	AcdTilePmt(1, 221),  AcdTilePmt(1,  20),
	AcdTilePmt(1,  10),  AcdTilePmt(1,   0),
	AcdTilePmt(1, 200),  AcdTilePmt(1, 1006),
	AcdTilePmt(1, 210),  AcdTilePmt(1, 220),
	AcdTilePmt(1, 1003),  AcdTilePmt(1, 230)},
       {AcdTilePmt(0, 222),  AcdTilePmt(0, 602), // 4
	AcdTilePmt(0, 212),  AcdTilePmt(0, 202),
	AcdTilePmt(0,  12),  AcdTilePmt(0,  13),
	AcdTilePmt(0,   3),  AcdTilePmt(0, 203),
	AcdTilePmt(0, 213),  AcdTilePmt(0, 223),
	AcdTilePmt(0,  14),  AcdTilePmt(0,   4),
	AcdTilePmt(0, 204),  AcdTilePmt(0, 1007),
	AcdTilePmt(0, 214),  AcdTilePmt(0, 224),
	AcdTilePmt(0, 1010), AcdTilePmt(0, 230)},
       {AcdTilePmt(1, 310),  AcdTilePmt(1, 300), // 5
	AcdTilePmt(1, 224),  AcdTilePmt(1, 214),
	AcdTilePmt(1, 603),  AcdTilePmt(1, 204),
	AcdTilePmt(1,   4),  AcdTilePmt(1,  14),
	AcdTilePmt(1, 223),  AcdTilePmt(1, 213),
	AcdTilePmt(1, 203),  AcdTilePmt(1,   3),
	AcdTilePmt(1,  13),  AcdTilePmt(1,  12),
	AcdTilePmt(1, 202),  AcdTilePmt(1, 212),
	AcdTilePmt(1, 1007),  AcdTilePmt(1, 222)},
       {AcdTilePmt(0, 324),  AcdTilePmt(0, 1009), // 6
	AcdTilePmt(0, 501),  AcdTilePmt(0, 1008),
	AcdTilePmt(0, 323),  AcdTilePmt(0, 313),
	AcdTilePmt(0, 303),  AcdTilePmt(0, 322),
	AcdTilePmt(0, 312),  AcdTilePmt(0, 302),
	AcdTilePmt(0, 301),  AcdTilePmt(0, 311),
	AcdTilePmt(0, 321),  AcdTilePmt(0, 500),
	AcdTilePmt(0, 300),  AcdTilePmt(0, 310),
	AcdTilePmt(0, 320),  AcdTilePmt(0, 330)},
       {AcdTilePmt(1, 330),  AcdTilePmt(1, 324), // 7
	AcdTilePmt(1, 314),  AcdTilePmt(1, 304),
	AcdTilePmt(1, 503),  AcdTilePmt(1, 323),
	AcdTilePmt(1, 313),  AcdTilePmt(1, 303),
	AcdTilePmt(1, 302),  AcdTilePmt(1, 312),
	AcdTilePmt(1, 322),  AcdTilePmt(1, 301),
	AcdTilePmt(1, 311),  AcdTilePmt(1, 321),
	AcdTilePmt(1, 1008),  AcdTilePmt(1, 502),
	AcdTilePmt(1, 1009),  AcdTilePmt(1, 320)},
       {AcdTilePmt(0, 314),  AcdTilePmt(0, 304), // 8
	AcdTilePmt(0, 1000),  AcdTilePmt(0, 424),
	AcdTilePmt(0, 603),  AcdTilePmt(0, 414),
	AcdTilePmt(0, 404),  AcdTilePmt(0,  44),
	AcdTilePmt(0,  34),  AcdTilePmt(0,  24),
	AcdTilePmt(0, 423),  AcdTilePmt(0, 413),
	AcdTilePmt(0, 403),  AcdTilePmt(0,  43),
	AcdTilePmt(0,  33),  AcdTilePmt(0,  23),
	AcdTilePmt(0, 1001),  AcdTilePmt(0,  42)},
       {AcdTilePmt(1,  42),  AcdTilePmt(1, 602), // 9
	AcdTilePmt(1,  23),  AcdTilePmt(1,  33),
	AcdTilePmt(1,  43),  AcdTilePmt(1, 403),
	AcdTilePmt(1, 413),  AcdTilePmt(1, 423),
	AcdTilePmt(1,  24),  AcdTilePmt(1,  34),
	AcdTilePmt(1,  44),  AcdTilePmt(1, 404),
	AcdTilePmt(1, 414),  AcdTilePmt(1, 1010),
	AcdTilePmt(1, 424),  AcdTilePmt(1, 1000),
	AcdTilePmt(1, 1001),  AcdTilePmt(1, 430)},
       {AcdTilePmt(0, 422),  AcdTilePmt(0, 601), // 10
	AcdTilePmt(0, 412),  AcdTilePmt(0, 402),
	AcdTilePmt(0,  32),  AcdTilePmt(0,  31),
	AcdTilePmt(0,  41),  AcdTilePmt(0, 401),
	AcdTilePmt(0, 411),  AcdTilePmt(0, 421),
	AcdTilePmt(0,  30),  AcdTilePmt(0,  40),
	AcdTilePmt(0, 400),  AcdTilePmt(0, 1002),
	AcdTilePmt(0, 410),  AcdTilePmt(0, 420),
	AcdTilePmt(0, 1003),  AcdTilePmt(0, 430)},
       {AcdTilePmt(1, 114),  AcdTilePmt(1, 104), // 11
	AcdTilePmt(1, 420),  AcdTilePmt(1, 410),
	AcdTilePmt(1, 600),  AcdTilePmt(1, 400),
	AcdTilePmt(1,  40),  AcdTilePmt(1,  30),
	AcdTilePmt(1, 421),  AcdTilePmt(1, 411),
	AcdTilePmt(1, 401),  AcdTilePmt(1,  41),
	AcdTilePmt(1,  31),  AcdTilePmt(1,  32),
	AcdTilePmt(1, 402),  AcdTilePmt(1, 412),
	AcdTilePmt(1, 1002),  AcdTilePmt(1, 422)}};

    return &(a[0][0]);
  };


  /*!
   * \brief Private function for initializing the ACD gem index map.
   *
   *  This comment is copied from the function above:
   *
   *  This function really should be a private method in ACDmap but Visual Studio
   *  gives an Internal Compiler Error when compiling LDFdump.cpp.
   *
   */
  inline unsigned* _initGemTbl() {
    static unsigned a[128] = {200,201,202,203,204,   // X - minus side 0-15
			      210,211,212,213,214,
			      220,221,222,223,224,
			      230,
			      400,401,402,403,404, // X - plus side  16-31
			      410,411,412,413,414,
			      420,421,422,423,424,
			      430,
			      100,101,102,103,104, // Y - minus side 32-47
			      110,111,112,113,114,
			      120,121,122,123,124,
			      130,
			      300,301,302,303,304, // Y - plus side 48-63
			      310,311,312,313,314,
			      320,321,322,323,324,
			      330,
			      000,001,002,003,004, // Top 64-88
			      010,011,012,013,014,
			      020,021,022,023,024,
			      030,031,032,033,034,
			      040,041,042,043,044,
			      0xFFFF,0xFFFF,0xFFFF,0xFFFF,0xFFFF, // 89-95 Must Be Zero
			      0xFFFF,0xFFFF,
			      500,501,502,503,       // Ribbons 96-103 
			      600,601,602,603,
			      0xFFFF,0xFFFF,0xFFFF,0xFFFF,       // 104 - 111 Must Be Zero
			      0xFFFF,0xFFFF,0xFFFF,0xFFFF,
			      1000,1001,1002,1003,1004,1005, // Not Assigned 112-122 
			      1006,1007,1008,1009,1010,  
			      0xFFFF,0xFFFF,0xFFFF,0xFFFF,0xFFFF};       // 123-127 Must Be Zero    
    return a;
  };
  
  const std::map<AcdTilePmt,AcdGarcGafe> _initMap() {
    static std::map<AcdTilePmt,AcdGarcGafe> theMap;
    if ( theMap.size() == 0 ) {
      for ( unsigned garc(0); garc < 12; garc++ ) {
	for ( unsigned gafe(0); gafe < 18; gafe++ ) {
	  const AcdTilePmt* tilePmt = toTilePmt(garc,gafe);
	  theMap[*tilePmt] = AcdGarcGafe(garc,gafe);
	}
      }
    }
    return theMap;
  };
  
  const std::map<unsigned,unsigned> _initGemMap() {   
    static std::map<unsigned,unsigned> theMap;
    if ( theMap.size() == 0 ) {
      for ( unsigned idx(0); idx < 128; idx++ ) {
	unsigned tile = tileFromGemIndex(idx);
	if ( tile != 0xFFFF ) {
	  theMap[tile] = idx;
	}
      }
    }
    return theMap;
  };

}

#endif // __CINT__


#endif
