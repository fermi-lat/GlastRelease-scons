# Microsoft Developer Studio Generated NMAKE File, Based on geometrylib.dsp
!IF "$(CFG)" == ""
CFG=geometrylib - Win32 Debug
!MESSAGE No configuration specified. Defaulting to geometrylib - Win32 Debug.
!ENDIF 

!IF "$(CFG)" != "geometrylib - Win32 Release" && "$(CFG)" != "geometrylib - Win32 Debug"
!MESSAGE Invalid configuration "$(CFG)" specified.
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "geometrylib.mak" CFG="geometrylib - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "geometrylib - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "geometrylib - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 
!ERROR An invalid configuration is specified.
!ENDIF 

!IF "$(OS)" == "Windows_NT"
NULL=
!ELSE 
NULL=nul
!ENDIF 

CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "geometrylib - Win32 Release"

OUTDIR=.\../../lib/Release
INTDIR=.\../../tmp/Release/geometry
# Begin Custom Macros
OutDir=.\../../lib/Release
# End Custom Macros

ALL : "$(OUTDIR)\geometry.lib"


CLEAN :
	-@erase "$(INTDIR)\Box.obj"
	-@erase "$(INTDIR)\Cone.obj"
	-@erase "$(INTDIR)\Conic.obj"
	-@erase "$(INTDIR)\CoordTransform.obj"
	-@erase "$(INTDIR)\Cylinder.obj"
	-@erase "$(INTDIR)\GeomObject.obj"
	-@erase "$(INTDIR)\Helix.obj"
	-@erase "$(INTDIR)\Hexahedron.obj"
	-@erase "$(INTDIR)\Hype.obj"
	-@erase "$(INTDIR)\Hyperbolic.obj"
	-@erase "$(INTDIR)\Intersection.obj"
	-@erase "$(INTDIR)\Plane.obj"
	-@erase "$(INTDIR)\Point.obj"
	-@erase "$(INTDIR)\Ray.obj"
	-@erase "$(INTDIR)\Shape.obj"
	-@erase "$(INTDIR)\Sphe.obj"
	-@erase "$(INTDIR)\Sphere.obj"
	-@erase "$(INTDIR)\Surface.obj"
	-@erase "$(INTDIR)\Track.obj"
	-@erase "$(INTDIR)\Tube.obj"
	-@erase "$(INTDIR)\vc60.idb"
	-@erase "$(INTDIR)\Vector.obj"
	-@erase "$(INTDIR)\Volume.obj"
	-@erase "$(INTDIR)\Wedge.obj"
	-@erase "$(OUTDIR)\geometry.lib"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

"$(INTDIR)" :
    if not exist "$(INTDIR)/$(NULL)" mkdir "$(INTDIR)"

CPP_PROJ=/nologo /G6 /ML /W3 /GR /GX /O2 /I "../" /I "../src" /D "WIN32" /D "_MBCS" /D "NDEBUG" /D "_LIB" /Fp"$(INTDIR)\geometrylib.pch" /YX /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /FI../mgr/nt_config.h /I../ @fincludes /c 
BSC32=bscmake.exe
BSC32_FLAGS=/nologo /o"$(OUTDIR)\geometrylib.bsc" 
BSC32_SBRS= \
	
LIB32=link.exe -lib
LIB32_FLAGS=/nologo /out:"$(OUTDIR)\geometry.lib" 
LIB32_OBJS= \
	"$(INTDIR)\Box.obj" \
	"$(INTDIR)\Cone.obj" \
	"$(INTDIR)\Conic.obj" \
	"$(INTDIR)\CoordTransform.obj" \
	"$(INTDIR)\Cylinder.obj" \
	"$(INTDIR)\GeomObject.obj" \
	"$(INTDIR)\Helix.obj" \
	"$(INTDIR)\Hexahedron.obj" \
	"$(INTDIR)\Hype.obj" \
	"$(INTDIR)\Hyperbolic.obj" \
	"$(INTDIR)\Intersection.obj" \
	"$(INTDIR)\Plane.obj" \
	"$(INTDIR)\Point.obj" \
	"$(INTDIR)\Ray.obj" \
	"$(INTDIR)\Shape.obj" \
	"$(INTDIR)\Sphe.obj" \
	"$(INTDIR)\Sphere.obj" \
	"$(INTDIR)\Surface.obj" \
	"$(INTDIR)\Track.obj" \
	"$(INTDIR)\Tube.obj" \
	"$(INTDIR)\Vector.obj" \
	"$(INTDIR)\Volume.obj" \
	"$(INTDIR)\Wedge.obj"

"$(OUTDIR)\geometry.lib" : "$(OUTDIR)" $(DEF_FILE) $(LIB32_OBJS)
    $(LIB32) @<<
  $(LIB32_FLAGS) $(DEF_FLAGS) $(LIB32_OBJS)
<<

!ELSEIF  "$(CFG)" == "geometrylib - Win32 Debug"

OUTDIR=.\../../lib/Debug
INTDIR=.\../../tmp/Debug/geometry
# Begin Custom Macros
OutDir=.\../../lib/Debug
# End Custom Macros

ALL : "$(OUTDIR)\geometry.lib"


CLEAN :
	-@erase "$(INTDIR)\Box.obj"
	-@erase "$(INTDIR)\Cone.obj"
	-@erase "$(INTDIR)\Conic.obj"
	-@erase "$(INTDIR)\CoordTransform.obj"
	-@erase "$(INTDIR)\Cylinder.obj"
	-@erase "$(INTDIR)\GeomObject.obj"
	-@erase "$(INTDIR)\Helix.obj"
	-@erase "$(INTDIR)\Hexahedron.obj"
	-@erase "$(INTDIR)\Hype.obj"
	-@erase "$(INTDIR)\Hyperbolic.obj"
	-@erase "$(INTDIR)\Intersection.obj"
	-@erase "$(INTDIR)\Plane.obj"
	-@erase "$(INTDIR)\Point.obj"
	-@erase "$(INTDIR)\Ray.obj"
	-@erase "$(INTDIR)\Shape.obj"
	-@erase "$(INTDIR)\Sphe.obj"
	-@erase "$(INTDIR)\Sphere.obj"
	-@erase "$(INTDIR)\Surface.obj"
	-@erase "$(INTDIR)\Track.obj"
	-@erase "$(INTDIR)\Tube.obj"
	-@erase "$(INTDIR)\vc60.idb"
	-@erase "$(INTDIR)\vc60.pdb"
	-@erase "$(INTDIR)\Vector.obj"
	-@erase "$(INTDIR)\Volume.obj"
	-@erase "$(INTDIR)\Wedge.obj"
	-@erase "$(OUTDIR)\geometry.lib"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

"$(INTDIR)" :
    if not exist "$(INTDIR)/$(NULL)" mkdir "$(INTDIR)"

CPP_PROJ=/nologo /G6 /MLd /W3 /Gm /GR /GX /ZI /Od /I "../" /I "../src" /D "WIN32" /D "_MBCS" /D "_DEBUG" /D "_LIB" /Fp"$(INTDIR)\geometrylib.pch" /YX /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /FI../mgr/nt_config.h /I../ /GZ @fincludes /c 
BSC32=bscmake.exe
BSC32_FLAGS=/nologo /o"$(OUTDIR)\geometrylib.bsc" 
BSC32_SBRS= \
	
LIB32=link.exe -lib
LIB32_FLAGS=/nologo /out:"$(OUTDIR)\geometry.lib" 
LIB32_OBJS= \
	"$(INTDIR)\Box.obj" \
	"$(INTDIR)\Cone.obj" \
	"$(INTDIR)\Conic.obj" \
	"$(INTDIR)\CoordTransform.obj" \
	"$(INTDIR)\Cylinder.obj" \
	"$(INTDIR)\GeomObject.obj" \
	"$(INTDIR)\Helix.obj" \
	"$(INTDIR)\Hexahedron.obj" \
	"$(INTDIR)\Hype.obj" \
	"$(INTDIR)\Hyperbolic.obj" \
	"$(INTDIR)\Intersection.obj" \
	"$(INTDIR)\Plane.obj" \
	"$(INTDIR)\Point.obj" \
	"$(INTDIR)\Ray.obj" \
	"$(INTDIR)\Shape.obj" \
	"$(INTDIR)\Sphe.obj" \
	"$(INTDIR)\Sphere.obj" \
	"$(INTDIR)\Surface.obj" \
	"$(INTDIR)\Track.obj" \
	"$(INTDIR)\Tube.obj" \
	"$(INTDIR)\Vector.obj" \
	"$(INTDIR)\Volume.obj" \
	"$(INTDIR)\Wedge.obj"

"$(OUTDIR)\geometry.lib" : "$(OUTDIR)" $(DEF_FILE) $(LIB32_OBJS)
    $(LIB32) @<<
  $(LIB32_FLAGS) $(DEF_FLAGS) $(LIB32_OBJS)
<<

!ENDIF 

.c{$(INTDIR)}.obj::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cpp{$(INTDIR)}.obj::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cxx{$(INTDIR)}.obj::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.c{$(INTDIR)}.sbr::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cpp{$(INTDIR)}.sbr::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cxx{$(INTDIR)}.sbr::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<


!IF "$(NO_EXTERNAL_DEPS)" != "1"
!IF EXISTS("geometrylib.dep")
!INCLUDE "geometrylib.dep"
!ELSE 
!MESSAGE Warning: cannot find "geometrylib.dep"
!ENDIF 
!ENDIF 


!IF "$(CFG)" == "geometrylib - Win32 Release" || "$(CFG)" == "geometrylib - Win32 Debug"
SOURCE=../src/Box.cxx

"$(INTDIR)\Box.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=../src/Cone.cxx

"$(INTDIR)\Cone.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=../src/Conic.cxx

"$(INTDIR)\Conic.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=../src/CoordTransform.cxx

"$(INTDIR)\CoordTransform.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=../src/Cylinder.cxx

"$(INTDIR)\Cylinder.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=../src/GeomObject.cxx

"$(INTDIR)\GeomObject.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=../src/Helix.cxx

"$(INTDIR)\Helix.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=../src/Hexahedron.cxx

"$(INTDIR)\Hexahedron.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=../src/Hype.cxx

"$(INTDIR)\Hype.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=../src/Hyperbolic.cxx

"$(INTDIR)\Hyperbolic.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=../src/Intersection.cxx

"$(INTDIR)\Intersection.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=../src/Plane.cxx

"$(INTDIR)\Plane.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=../src/Point.cxx

"$(INTDIR)\Point.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=../src/Ray.cxx

"$(INTDIR)\Ray.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=../src/Shape.cxx

"$(INTDIR)\Shape.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=../src/Sphe.cxx

"$(INTDIR)\Sphe.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=../src/Sphere.cxx

"$(INTDIR)\Sphere.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=../src/Surface.cxx

"$(INTDIR)\Surface.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=../src/Track.cxx

"$(INTDIR)\Track.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=../src/Tube.cxx

"$(INTDIR)\Tube.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=../src/Vector.cxx

"$(INTDIR)\Vector.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=../src/Volume.cxx

"$(INTDIR)\Volume.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=../src/Wedge.cxx

"$(INTDIR)\Wedge.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)



!ENDIF 

