# Microsoft Developer Studio Project File - Name="geomreplib" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=geomreplib - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "geomreplib.mak".
!MESSAGE 
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

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "geomreplib - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "../Win32"
# PROP BASE Intermediate_Dir "../Win32"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "../../msdev/lib"
# PROP Intermediate_Dir "../../msdev/tmp/geomrep"
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
# ADD LIB32 /nologo /out:"../../msdev/lib/geomrep.lib"

!ELSEIF  "$(CFG)" == "geomreplib - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "../Win32Debug"
# PROP BASE Intermediate_Dir "../Win32Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "../../msdev/Debug/lib"
# PROP Intermediate_Dir "../../msdev/Debug/tmp/geomrep"
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
# ADD LIB32 /nologo /out:"../../msdev/Debug/lib/geomrep.lib"

!ENDIF 

# Begin Target

# Name "geomreplib - Win32 Release"
# Name "geomreplib - Win32 Debug"
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

SOURCE=../src/ArcRep.cxx
# End Source File
# Begin Source File

SOURCE=../geomrep/ArcRep.h
# End Source File
# Begin Source File

SOURCE=../src/BoxRep.cxx
# End Source File
# Begin Source File

SOURCE=../geomrep/BoxRep.h
# End Source File
# Begin Source File

SOURCE=../src/ConeRep.cxx
# End Source File
# Begin Source File

SOURCE=../geomrep/ConeRep.h
# End Source File
# Begin Source File

SOURCE=../src/ConsRep.cxx
# End Source File
# Begin Source File

SOURCE=../geomrep/ConsRep.h
# End Source File
# Begin Source File

SOURCE=../src/HypeRep.cxx
# End Source File
# Begin Source File

SOURCE=../geomrep/HypeRep.h
# End Source File
# Begin Source File

SOURCE=../src/PointRep.cxx
# End Source File
# Begin Source File

SOURCE=../geomrep/PointRep.h
# End Source File
# Begin Source File

SOURCE=../src/SpheRep.cxx
# End Source File
# Begin Source File

SOURCE=../geomrep/SpheRep.h
# End Source File
# Begin Source File

SOURCE=../src/TrackRep.cxx
# End Source File
# Begin Source File

SOURCE=../geomrep/TrackRep.h
# End Source File
# Begin Source File

SOURCE=../src/TubeRep.cxx
# End Source File
# Begin Source File

SOURCE=../geomrep/TubeRep.h
# End Source File
# Begin Source File

SOURCE=../src/TubsRep.cxx
# End Source File
# Begin Source File

SOURCE=../geomrep/TubsRep.h
# End Source File
# End Group
# End Target
# End Project
