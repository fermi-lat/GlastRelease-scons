# Microsoft Developer Studio Project File - Name="geometrylib" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=geometrylib - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "geometrylib.mak".
!MESSAGE 
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

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "geometrylib - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "../Win32"
# PROP BASE Intermediate_Dir "../Win32"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "../../lib/Release"
# PROP Intermediate_Dir "../../tmp/Release/geometry"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /MD /W3 /GR /GX /O2 /I "../" /I "../src" /D "WIN32" /D "_MBCS" /D "NDEBUG" /D "_LIB" /YX /FD /FI../mgr/nt_config.h @fincludes /c
# ADD CPP /nologo /G6 /W3 /GR /GX /O2 /I "../" /I "../src" /D "WIN32" /D "_MBCS" /D "NDEBUG" /D "_LIB" /YX /FD /FI../mgr/nt_config.h /I../ @fincludes /c
# ADD BASE RSC /l 0x40c /d "NDEBUG"
# ADD RSC /l 0x40c /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo /out:"../../lib/Release/geometry.lib"

!ELSEIF  "$(CFG)" == "geometrylib - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "../Win32Debug"
# PROP BASE Intermediate_Dir "../Win32Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "../../lib/Debug"
# PROP Intermediate_Dir "../../tmp/Debug/geometry"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /MDd /W3 /Gm /GR /GX /ZI /Od /I "../" /I "../src" /D "WIN32" /D "_MBCS" /D "_DEBUG" /D "_LIB" /YX /FD /FI../mgr/nt_config.h /GZ @fincludes /c
# ADD CPP /nologo /G6 /W3 /Gm /GR /GX /ZI /Od /I "../" /I "../src" /D "WIN32" /D "_MBCS" /D "_DEBUG" /D "_LIB" /YX /FD /FI../mgr/nt_config.h /I../ /GZ @fincludes /c
# ADD BASE RSC /l 0x40c /d "_DEBUG"
# ADD RSC /l 0x40c /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo /out:"../../lib/Debug/geometry.lib"

!ENDIF 

# Begin Target

# Name "geometrylib - Win32 Release"
# Name "geometrylib - Win32 Debug"
# Begin Group "config"

# PROP Default_Filter ""
# Begin Source File

SOURCE=fincludes
# End Source File
# Begin Source File

SOURCE=..\src\Makefile.am
# End Source File
# End Group
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE=../src/Box.cxx
# End Source File
# Begin Source File

SOURCE=../geometry/Box.h
# End Source File
# Begin Source File

SOURCE=../src/Cone.cxx
# End Source File
# Begin Source File

SOURCE=../geometry/Cone.h
# End Source File
# Begin Source File

SOURCE=../src/Conic.cxx
# End Source File
# Begin Source File

SOURCE=../geometry/Conic.h
# End Source File
# Begin Source File

SOURCE=../geometry/Cons.h
# End Source File
# Begin Source File

SOURCE=../geometry/CoordSystem.h
# End Source File
# Begin Source File

SOURCE=../src/CoordTransform.cxx
# End Source File
# Begin Source File

SOURCE=../geometry/CoordTransform.h
# End Source File
# Begin Source File

SOURCE=../src/Cylinder.cxx
# End Source File
# Begin Source File

SOURCE=../geometry/Cylinder.h
# End Source File
# Begin Source File

SOURCE=../geometry/Digitizer.h
# End Source File
# Begin Source File

SOURCE=../src/GeomObject.cxx
# End Source File
# Begin Source File

SOURCE=../geometry/GeomObject.h
# End Source File
# Begin Source File

SOURCE=../src/Helix.cxx
# End Source File
# Begin Source File

SOURCE=../geometry/Helix.h
# End Source File
# Begin Source File

SOURCE=../src/Hexahedron.cxx
# End Source File
# Begin Source File

SOURCE=../geometry/Hexahedron.h
# End Source File
# Begin Source File

SOURCE=../src/Hype.cxx
# End Source File
# Begin Source File

SOURCE=../geometry/Hype.h
# End Source File
# Begin Source File

SOURCE=../src/Hyperbolic.cxx
# End Source File
# Begin Source File

SOURCE=../geometry/Hyperbolic.h
# End Source File
# Begin Source File

SOURCE=../src/Intersection.cxx
# End Source File
# Begin Source File

SOURCE=../geometry/Intersection.h
# End Source File
# Begin Source File

SOURCE=../src/Plane.cxx
# End Source File
# Begin Source File

SOURCE=../geometry/Plane.h
# End Source File
# Begin Source File

SOURCE=../src/Point.cxx
# End Source File
# Begin Source File

SOURCE=../geometry/Point.h
# End Source File
# Begin Source File

SOURCE=../src/Ray.cxx
# End Source File
# Begin Source File

SOURCE=../geometry/Ray.h
# End Source File
# Begin Source File

SOURCE=../src/Shape.cxx
# End Source File
# Begin Source File

SOURCE=../geometry/Shape.h
# End Source File
# Begin Source File

SOURCE=../src/Sphe.cxx
# End Source File
# Begin Source File

SOURCE=../geometry/Sphe.h
# End Source File
# Begin Source File

SOURCE=../src/Sphere.cxx
# End Source File
# Begin Source File

SOURCE=../geometry/Sphere.h
# End Source File
# Begin Source File

SOURCE=../src/Surface.cxx
# End Source File
# Begin Source File

SOURCE=../geometry/Surface.h
# End Source File
# Begin Source File

SOURCE=../src/Track.cxx
# End Source File
# Begin Source File

SOURCE=../geometry/Track.h
# End Source File
# Begin Source File

SOURCE=../geometry/Trap.h
# End Source File
# Begin Source File

SOURCE=../geometry/Trd1.h
# End Source File
# Begin Source File

SOURCE=../geometry/Trd2.h
# End Source File
# Begin Source File

SOURCE=../src/Tube.cxx
# End Source File
# Begin Source File

SOURCE=../geometry/Tube.h
# End Source File
# Begin Source File

SOURCE=../geometry/Tubs.h
# End Source File
# Begin Source File

SOURCE=../src/Vector.cxx
# End Source File
# Begin Source File

SOURCE=../geometry/Vector.h
# End Source File
# Begin Source File

SOURCE=../src/Volume.cxx
# End Source File
# Begin Source File

SOURCE=../geometry/Volume.h
# End Source File
# Begin Source File

SOURCE=../src/Wedge.cxx
# End Source File
# Begin Source File

SOURCE=../geometry/Wedge.h
# End Source File
# End Group
# End Target
# End Project
