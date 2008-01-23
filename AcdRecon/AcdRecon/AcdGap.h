#ifndef ACDGAP_H
#define ACDGAP_H

namespace AcdRecon {

  typedef enum { None=0,             
		 X_RibbonSide=1,     
		 Y_RibbonSide=2, 
		 Y_RibbonTop=3,
		 SideCornerEdge=4,
		 TopCornerEdge=5,
		 TileHole=6,
		 CornerRay=7,
		 NumGapTypes } AcdGapType;

  typedef enum { Top=0,
		 MinusX,
		 MimusY,
		 PlusX,
		 PluxY,
		 Bottom,
		 NumFaces } AcdFace;

}


#endif
