# Microsoft Developer Studio Generated NMAKE File, Based on geomreplib.dsp
!IF "$(CFG)" == ""
CFG=geomreplib - Win32 Debug
!MESSAGE No configuration specified. Defaulting to geomreplib - Win32 Debug.
!ENDIF 

!IF "$(CFG)" != "geomreplib - Win32 Release" && "$(CFG)" != "geomreplib - Win32 Debug"
!MESSAGE Invalid configuration "$(CFG)" specified.
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "geomreplib.mak" CFG="geomreplib - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "geomreplib - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "geomreplib - Win32 Debug" (based on "Win32 (x86) Static Library")
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

!IF  "$(CFG)" == "geomreplib - Win32 Release"

OUTDIR=.\../../lib/Release
INTDIR=.\../../tmp/Release/geomrep
# Begin Custom Macros
OutDir=.\../../lib/Release
# End Custom Macros

ALL : "$(OUTDIR)\geomrep.lib"


CLEAN :
	-@erase "$(INTDIR)\ArcRep.obj"
	-@erase "$(INTDIR)\BoxRep.obj"
	-@erase "$(INTDIR)\ConeRep.obj"
	-@erase "$(INTDIR)\ConsRep.obj"
	-@erase "$(INTDIR)\HypeRep.obj"
	-@erase "$(INTDIR)\PointRep.obj"
	-@erase "$(INTDIR)\SpheRep.obj"
	-@erase "$(INTDIR)\TrackRep.obj"
	-@erase "$(INTDIR)\TubeRep.obj"
	-@erase "$(INTDIR)\TubsRep.obj"
	-@erase "$(INTDIR)\vc60.idb"
	-@erase "$(OUTDIR)\geomrep.lib"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

"$(INTDIR)" :
    if not exist "$(INTDIR)/$(NULL)" mkdir "$(INTDIR)"

CPP_PROJ=/nologo /G6 /ML /W3 /GR /GX /O2 /I "../" /I "../src" /D "WIN32" /D "_MBCS" /D "NDEBUG" /D "_LIB" /Fp"$(INTDIR)\geomreplib.pch" /YX /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /FI../mgr/nt_config.h /I../ @fincludes /c 
BSC32=bscmake.exe
BSC32_FLAGS=/nologo /o"$(OUTDIR)\geomreplib.bsc" 
BSC32_SBRS= \
	
LIB32=link.exe -lib
LIB32_FLAGS=/nologo /out:"$(OUTDIR)\geomrep.lib" 
LIB32_OBJS= \
	"$(INTDIR)\ArcRep.obj" \
	"$(INTDIR)\BoxRep.obj" \
	"$(INTDIR)\ConeRep.obj" \
	"$(INTDIR)\ConsRep.obj" \
	"$(INTDIR)\HypeRep.obj" \
	"$(INTDIR)\PointRep.obj" \
	"$(INTDIR)\SpheRep.obj" \
	"$(INTDIR)\TrackRep.obj" \
	"$(INTDIR)\TubeRep.obj" \
	"$(INTDIR)\TubsRep.obj"

"$(OUTDIR)\geomrep.lib" : "$(OUTDIR)" $(DEF_FILE) $(LIB32_OBJS)
    $(LIB32) @<<
  $(LIB32_FLAGS) $(DEF_FLAGS) $(LIB32_OBJS)
<<

!ELSEIF  "$(CFG)" == "geomreplib - Win32 Debug"

OUTDIR=.\../../lib/Debug
INTDIR=.\../../tmp/Debug/geomrep
# Begin Custom Macros
OutDir=.\../../lib/Debug
# End Custom Macros

ALL : "$(OUTDIR)\geomrep.lib"


CLEAN :
	-@erase "$(INTDIR)\ArcRep.obj"
	-@erase "$(INTDIR)\BoxRep.obj"
	-@erase "$(INTDIR)\ConeRep.obj"
	-@erase "$(INTDIR)\ConsRep.obj"
	-@erase "$(INTDIR)\HypeRep.obj"
	-@erase "$(INTDIR)\PointRep.obj"
	-@erase "$(INTDIR)\SpheRep.obj"
	-@erase "$(INTDIR)\TrackRep.obj"
	-@erase "$(INTDIR)\TubeRep.obj"
	-@erase "$(INTDIR)\TubsRep.obj"
	-@erase "$(INTDIR)\vc60.idb"
	-@erase "$(INTDIR)\vc60.pdb"
	-@erase "$(OUTDIR)\geomrep.lib"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

"$(INTDIR)" :
    if not exist "$(INTDIR)/$(NULL)" mkdir "$(INTDIR)"

CPP_PROJ=/nologo /G6 /MLd /W3 /Gm /GR /GX /ZI /Od /I "../" /I "../src" /D "WIN32" /D "_MBCS" /D "_DEBUG" /D "_LIB" /Fp"$(INTDIR)\geomreplib.pch" /YX /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /FI../mgr/nt_config.h /I../ /GZ @fincludes /c 
BSC32=bscmake.exe
BSC32_FLAGS=/nologo /o"$(OUTDIR)\geomreplib.bsc" 
BSC32_SBRS= \
	
LIB32=link.exe -lib
LIB32_FLAGS=/nologo /out:"$(OUTDIR)\geomrep.lib" 
LIB32_OBJS= \
	"$(INTDIR)\ArcRep.obj" \
	"$(INTDIR)\BoxRep.obj" \
	"$(INTDIR)\ConeRep.obj" \
	"$(INTDIR)\ConsRep.obj" \
	"$(INTDIR)\HypeRep.obj" \
	"$(INTDIR)\PointRep.obj" \
	"$(INTDIR)\SpheRep.obj" \
	"$(INTDIR)\TrackRep.obj" \
	"$(INTDIR)\TubeRep.obj" \
	"$(INTDIR)\TubsRep.obj"

"$(OUTDIR)\geomrep.lib" : "$(OUTDIR)" $(DEF_FILE) $(LIB32_OBJS)
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
!IF EXISTS("geomreplib.dep")
!INCLUDE "geomreplib.dep"
!ELSE 
!MESSAGE Warning: cannot find "geomreplib.dep"
!ENDIF 
!ENDIF 


!IF "$(CFG)" == "geomreplib - Win32 Release" || "$(CFG)" == "geomreplib - Win32 Debug"
SOURCE=../src/ArcRep.cxx

"$(INTDIR)\ArcRep.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=../src/BoxRep.cxx

"$(INTDIR)\BoxRep.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=../src/ConeRep.cxx

"$(INTDIR)\ConeRep.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=../src/ConsRep.cxx

"$(INTDIR)\ConsRep.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=../src/HypeRep.cxx

"$(INTDIR)\HypeRep.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=../src/PointRep.cxx

"$(INTDIR)\PointRep.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=../src/SpheRep.cxx

"$(INTDIR)\SpheRep.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=../src/TrackRep.cxx

"$(INTDIR)\TrackRep.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=../src/TubeRep.cxx

"$(INTDIR)\TubeRep.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=../src/TubsRep.cxx

"$(INTDIR)\TubsRep.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)



!ENDIF 

