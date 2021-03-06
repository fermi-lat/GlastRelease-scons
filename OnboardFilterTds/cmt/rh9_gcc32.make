CMT_tag=$(tag)
CMTROOT=/afs/slac.stanford.edu/g/glast/applications/CMT/v1r16p20040701
CMT_root=/afs/slac.stanford.edu/g/glast/applications/CMT/v1r16p20040701
CMTVERSION=v1r16p20040701
CMTrelease=15
cmt_hardware_query_command=uname -m
cmt_hardware=`$(cmt_hardware_query_command)`
cmt_system_version_query_command=${CMTROOT}/mgr/cmt_linux_version.sh | ${CMTROOT}/mgr/cmt_filter_version.sh
cmt_system_version=`$(cmt_system_version_query_command)`
cmt_compiler_version_query_command=${CMTROOT}/mgr/cmt_gcc_version.sh | ${CMTROOT}/mgr/cmt_filter_version.sh
cmt_compiler_version=`$(cmt_compiler_version_query_command)`
PATH=${ROOTSYS}/bin:/afs/slac.stanford.edu/g/glast/applications/CMT/v1r16p20040701/${CMTBIN}:/nfs/slac/g/svac/borgland/Valgrind-3.1.0/bin:/nfs/slac/g/svac/borgland/valkyrie-1.1.0/bin:/u/gl/heather/bin:/usr/local/bin:/usr/afsws/bin:/usr/afsws/etc:/bin:/usr/bin:/usr/etc:/usr/bin/X11:.:/u/gl/heather/bin:/afs/slac.stanford.edu/g/glast/ground/scripts:/nfs/slac/g/glast/users/glground/heather/mpatrol/mpatrol/build/unix:/nfs/farm/g/glast/u06/berrie/astroroot/bin:/nfs/farm/g/glast/u06/berrie/astroroot/lib
CLASSPATH=/afs/slac.stanford.edu/g/glast/applications/CMT/v1r16p20040701/java
debug_option=-g
cc=gcc
cdebugflags=$(debug_option)
pp_cflags=-Di586
ccomp=$(cc) -c $(includes) $(cdebugflags) $(cflags) $(pp_cflags)
clink=$(cc) $(clinkflags) $(cdebugflags)
ppcmd=-I
preproc=c++ -MD -c 
cpp=g++
cppdebugflags=$(debug_option)
cppflags=-pipe -ansi -W -Wall  -fPIC -shared -D_GNU_SOURCE -Dlinux -Dunix 
pp_cppflags=-D_GNU_SOURCE
cppcomp=$(cpp) -c $(includes) $(cppoptions) $(cppflags) $(pp_cppflags)
cpplinkflags=-Wl,-Bdynamic  $(linkdebug)
cpplink=$(cpp)   $(cpplinkflags)
for=g77
fflags=$(debug_option)
fcomp=$(for) -c $(fincludes) $(fflags) $(pp_fflags)
flink=$(for) $(flinkflags)
javacomp=javac -classpath $(src):$(CLASSPATH) 
javacopy=cp
jar=jar
X11_cflags=-I/usr/include
Xm_cflags=-I/usr/include
X_linkopts=-L/usr/X11R6/lib -lXm -lXt -lXext -lX11 -lm
lex=flex $(lexflags)
yaccflags= -l -d 
yacc=yacc $(yaccflags)
ar=ar r
ranlib=ranlib
make_shlib=${CMTROOT}/mgr/cmt_make_shlib_common.sh extract
shlibsuffix=so
shlibbuilder=g++ $(cmt_installarea_linkopts) 
shlibflags=-shared
symlink=/bin/ln -fs 
symunlink=/bin/rm -f 
build_library_links=$(cmtexe) build library_links -quiet -tag=$(tags)
remove_library_links=$(cmtexe) remove library_links -quiet -tag=$(tags)
cmtexe=${CMTROOT}/${CMTBIN}/cmt.exe
build_prototype=$(cmtexe) build prototype
build_dependencies=$(cmtexe) -quiet -tag=$(tags) build dependencies
build_triggers=$(cmtexe) build triggers
implied_library_prefix=-l
SHELL=/bin/sh
src=../src/
doc=../doc/
inc=../src/
mgr=../cmt/
application_suffix=.exe
library_prefix=lib
lock_command=chmod -R a-w ../*
unlock_command=chmod -R g+w ../*
MAKEFLAGS= --no-print-directory 
gmake_hosts=lx1 rsplus lxtest as7 dxplus ax7 hp2 aleph hp1 hpplus papou1-fe atlas
make_hosts=virgo-control1 rio0a vmpc38a
everywhere=hosts
install_command=cp 
uninstall_command=/bin/rm -f 
cmt_installarea_command=ln -s 
cmt_uninstallarea_command=/bin/rm -f 
cmt_install_area_command=$(cmt_installarea_command)
cmt_uninstall_area_command=$(cmt_uninstallarea_command)
cmt_install_action=$(CMTROOT)/mgr/cmt_install_action.sh
cmt_installdir_action=$(CMTROOT)/mgr/cmt_installdir_action.sh
cmt_uninstall_action=$(CMTROOT)/mgr/cmt_uninstall_action.sh
cmt_uninstalldir_action=$(CMTROOT)/mgr/cmt_uninstalldir_action.sh
mkdir=mkdir
cmt_installarea_prefix=InstallArea
CMT_PATH_remove_regexp=/[^/]*/
CMT_PATH_remove_share_regexp=/share/
NEWCMTCONFIG=i686-rhel36-gcc32
OnboardFilterTds_tag=$(tag)
ONBOARDFILTERTDSROOT=/a/sulky28/g.glast.u10/builds/GlastRelease/gr-v9r10-build/OnboardFilterTds/v0
OnboardFilterTds_root=/a/sulky28/g.glast.u10/builds/GlastRelease/gr-v9r10-build/OnboardFilterTds/v0
ONBOARDFILTERTDSVERSION=v0
OnboardFilterTds_cmtpath=/nfs/farm/g/glast/u10/builds/GlastRelease/gr-v9r10-build
OnboardFilterTds_offset=a/sulky28/g.glast.u10/builds/GlastRelease/gr-v9r10-build
OnboardFilterTds_project=Project1
obf_tag=$(tag)
OBFROOT=/nfs/farm/g/glast/u10/builds/GlastRelease/gr-v9r10-build/IExternal/obf/v0r3
obf_root=/nfs/farm/g/glast/u10/builds/GlastRelease/gr-v9r10-build/IExternal/obf/v0r3
OBFVERSION=v0r3
obf_cmtpath=/nfs/farm/g/glast/u10/builds/GlastRelease/gr-v9r10-build
obf_offset=IExternal
obf_project=Project1
obf_native_version=v0r0
obf_DIR=${GLAST_EXT}/obf/${obf_native_version}
FLIGHTCODELIBS=/afs/slac/g/glast/ground/GLAST_EXT/rh9_gcc32/obf/v0r0/lib
includeSym=-I
includes=$(includeSym)${obf_DIR}/src              $(includeSym)${obf_DIR}/src/CAL_DB                  $(includeSym)${obf_DIR}/src/CDM                     $(includeSym)${obf_DIR}/src/CDM/src                 $(includeSym)${obf_DIR}/src/CGB_DB                  $(includeSym)${obf_DIR}/src/CMX                     $(includeSym)${obf_DIR}/src/CMX/src                 $(includeSym)${obf_DIR}/src/CPG_DB                  $(includeSym)${obf_DIR}/src/CPP_DB                  $(includeSym)${obf_DIR}/src/EDS                     $(includeSym)${obf_DIR}/src/EDS_DB                  $(includeSym)${obf_DIR}/src/EFC                     $(includeSym)${obf_DIR}/src/EFC_DB                  $(includeSym)${obf_DIR}/src/EMP                     $(includeSym)${obf_DIR}/src/GEO_DB                  $(includeSym)${obf_DIR}/src/GFC_DB                  $(includeSym)${obf_DIR}/src/GGF_DB                  $(includeSym)${obf_DIR}/src/MDB                     $(includeSym)${obf_DIR}/src/MPP_DB                  $(includeSym)${obf_DIR}/src/PBI                     $(includeSym)${obf_DIR}/src/PBS                     $(includeSym)${obf_DIR}/src/RFC_DB  $(ppcmd)"$(ONBOARDFILTERTDSROOT)" $(use_includes)
obf_linkopts= -L$(obf_DIR)/lib -lflight_cdm -lflight_cmx -lflight_mdb -leds -lefc  
LD_LIBRARY_PATH=/afs/slac/g/glast/ground/GLAST_EXT/rh9_gcc32/cppunit/1.10.2/lib:/afs/slac/g/glast/ground/GLAST_EXT/rh9_gcc32/CLHEP/1.9.2.2/lib:/afs/slac/g/glast/ground/GLAST_EXT/rh9_gcc32/ROOT/v5.10.00/root/bin:/afs/slac/g/glast/ground/GLAST_EXT/rh9_gcc32/ROOT/v5.10.00/root/lib:/afs/slac/g/glast/ground/GLAST_EXT/rh9_gcc32/obf/v0r0/lib:/afs/slac/g/glast/ground/GLAST_EXT/rh9_gcc32/xerces/2.6.0/lib:/afs/slac/g/glast/ground/GLAST_EXT/rh9_gcc32/gaudi/v18r1/lib:/nfs/farm/g/glast/u10/builds/GlastRelease/gr-v9r10-build/lib/
GaudiInterface_tag=$(tag)
GAUDIINTERFACEROOT=/nfs/farm/g/glast/u10/builds/GlastRelease/gr-v9r10-build/IExternal/GaudiInterface/v0r181p5
GaudiInterface_root=/nfs/farm/g/glast/u10/builds/GlastRelease/gr-v9r10-build/IExternal/GaudiInterface/v0r181p5
GAUDIINTERFACEVERSION=v0r181p5
GaudiInterface_cmtpath=/nfs/farm/g/glast/u10/builds/GlastRelease/gr-v9r10-build
GaudiInterface_offset=IExternal
GaudiInterface_project=Project1
ExternalLibs_tag=$(tag)
EXTERNALLIBSROOT=/nfs/farm/g/glast/u10/builds/GlastRelease/gr-v9r10-build/IExternal/ExternalLibs/v5r0
ExternalLibs_root=/nfs/farm/g/glast/u10/builds/GlastRelease/gr-v9r10-build/IExternal/ExternalLibs/v5r0
EXTERNALLIBSVERSION=v5r0
ExternalLibs_cmtpath=/nfs/farm/g/glast/u10/builds/GlastRelease/gr-v9r10-build
ExternalLibs_offset=IExternal
ExternalLibs_project=Project1
EXTPACK_DIR=${GLAST_EXT}
GaudiInterface_native_version=v18r1
GaudiInterface_DIR=${EXTPACK_DIR}/gaudi/${GaudiInterface_native_version}
gaudi_root=${GaudiInterface_DIR}
ROOT_tag=$(tag)
ROOTROOT=/nfs/farm/g/glast/u10/builds/GlastRelease/gr-v9r10-build/IExternal/ROOT/v5r10
ROOT_root=/nfs/farm/g/glast/u10/builds/GlastRelease/gr-v9r10-build/IExternal/ROOT/v5r10
ROOTVERSION=v5r10
ROOT_cmtpath=/nfs/farm/g/glast/u10/builds/GlastRelease/gr-v9r10-build
ROOT_offset=IExternal
ROOT_project=Project1
ROOT_DIR=${GLAST_EXT}/ROOT
ROOT_native_version=v5.10.00
ROOT_PATH=${ROOT_DIR}/$(ROOT_native_version)/root
ROOTSYS=/afs/slac/g/glast/ground/GLAST_EXT/rh9_gcc32/ROOT/v5.10.00/root
dict=../dict/
rootcint=rootcint
ROOT_libs=-L$(ROOT_PATH)/lib -lCore -lCint -lTree -lMatrix -lPhysics -lpthread -lm -ldl -rdynamic
ROOT_GUI_libs=-L$(ROOT_PATH)/lib -lHist -lGraf -lGraf3d -lGpad -lRint -lPostscript -lTreePlayer 
ROOT_linkopts=$(ROOT_libs)
ROOT_cppflagsEx=$(ppcmd) "$(ROOT_PATH)/include" -DUSE_ROOT
ROOT_cppflags=-fpermissive
CLHEP_tag=$(tag)
CLHEPROOT=/nfs/farm/g/glast/u10/builds/GlastRelease/gr-v9r10-build/IExternal/CLHEP/v3r0
CLHEP_root=/nfs/farm/g/glast/u10/builds/GlastRelease/gr-v9r10-build/IExternal/CLHEP/v3r0
CLHEPVERSION=v3r0
CLHEP_cmtpath=/nfs/farm/g/glast/u10/builds/GlastRelease/gr-v9r10-build
CLHEP_offset=IExternal
CLHEP_project=Project1
CLHEP_native_version=1.9.2.2
CLHEP_DIR=$(GLAST_EXT)/CLHEP
CLHEPBASE=${CLHEP_DIR}/$(CLHEP_native_version)
CLHEP_linkopts=-L$(CLHEPBASE)/lib -lCLHEP
cppunit_tag=$(tag)
CPPUNITROOT=/nfs/farm/g/glast/u10/builds/GlastRelease/gr-v9r10-build/IExternal/cppunit/v2r0
cppunit_root=/nfs/farm/g/glast/u10/builds/GlastRelease/gr-v9r10-build/IExternal/cppunit/v2r0
CPPUNITVERSION=v2r0
cppunit_cmtpath=/nfs/farm/g/glast/u10/builds/GlastRelease/gr-v9r10-build
cppunit_offset=IExternal
cppunit_project=Project1
cppunit_native_version=1.10.2
cppunit_DIR=${GLAST_EXT}/cppunit/$(cppunit_native_version)
cppunit_linkopts=-L ${cppunit_DIR}/lib/ -lcppunit -ldl 
XMLEXT_tag=$(tag)
XMLEXTROOT=/nfs/farm/g/glast/u10/builds/GlastRelease/gr-v9r10-build/IExternal/XMLEXT/v5r260p0
XMLEXT_root=/nfs/farm/g/glast/u10/builds/GlastRelease/gr-v9r10-build/IExternal/XMLEXT/v5r260p0
XMLEXTVERSION=v5r260p0
XMLEXT_cmtpath=/nfs/farm/g/glast/u10/builds/GlastRelease/gr-v9r10-build
XMLEXT_offset=IExternal
XMLEXT_project=Project1
XMLEXT_native_version=2.6.0
XMLEXT_DIR=${EXTPACK_DIR}/xerces/$(XMLEXT_native_version)
XMLEXT_linkopts=-L$(XMLEXT_DIR)/lib/ -lxerces-c -lpthread
BOOST_DIR=${GLAST_EXT}/boost
BOOSTHOME=$(BOOST_DIR)/1.31.0
BOOST_linkopts= -L${GaudiInterface_DIR}/lib -lboost_filesystem-gcc 
GaudiKernel_linkopts=-L$(GaudiInterface_DIR)/lib -lGaudiKernel -llcg_SealBase -llcg_PluginManager -llcg_SealKernel -llcg_Reflection -llcg_ReflectionBuilder -llcg_Reflex
GaudiPiShr=/afs/slac/g/glast/ground/GLAST_EXT/rh9_gcc32/gaudi/v18r1/lib/libGaudiPi
GaudiSvcShr=/afs/slac/g/glast/ground/GLAST_EXT/rh9_gcc32/gaudi/v18r1/lib/libGaudiSvc
GaudiAudShr=/afs/slac/g/glast/ground/GLAST_EXT/rh9_gcc32/gaudi/v18r1/lib/libGaudiAud
GaudiAlgShr=/afs/slac/g/glast/ground/GLAST_EXT/rh9_gcc32/gaudi/v18r1/lib/libGaudiAlg
RootHistCnvShr=/afs/slac/g/glast/ground/GLAST_EXT/rh9_gcc32/gaudi/v18r1/lib/libRootHistCnv
GaudiInterface_linkopts=$(GaudiKernel_linkopts) $(BOOST_linkopts)
Event_tag=$(tag)
EVENTROOT=/nfs/farm/g/glast/u10/builds/GlastRelease/gr-v9r10-build/Event/v11r13
Event_root=/nfs/farm/g/glast/u10/builds/GlastRelease/gr-v9r10-build/Event/v11r13
EVENTVERSION=v11r13
Event_cmtpath=/nfs/farm/g/glast/u10/builds/GlastRelease/gr-v9r10-build
Event_project=Project1
geometry_tag=$(tag)
GEOMETRYROOT=/nfs/farm/g/glast/u10/builds/GlastRelease/gr-v9r10-build/geometry/v3r2
geometry_root=/nfs/farm/g/glast/u10/builds/GlastRelease/gr-v9r10-build/geometry/v3r2
GEOMETRYVERSION=v3r2
geometry_cmtpath=/nfs/farm/g/glast/u10/builds/GlastRelease/gr-v9r10-build
geometry_project=Project1
GlastPolicy_tag=$(tag)
GLASTPOLICYROOT=/nfs/farm/g/glast/u10/builds/GlastRelease/gr-v9r10-build/GlastPolicy/v6r11
GlastPolicy_root=/nfs/farm/g/glast/u10/builds/GlastRelease/gr-v9r10-build/GlastPolicy/v6r11
GLASTPOLICYVERSION=v6r11
GlastPolicy_cmtpath=/nfs/farm/g/glast/u10/builds/GlastRelease/gr-v9r10-build
GlastPolicy_project=Project1
GlastPatternPolicy_tag=$(tag)
GLASTPATTERNPOLICYROOT=/nfs/farm/g/glast/u10/builds/GlastRelease/gr-v9r10-build/GlastPolicy/GlastPatternPolicy/v1r3p3
GlastPatternPolicy_root=/nfs/farm/g/glast/u10/builds/GlastRelease/gr-v9r10-build/GlastPolicy/GlastPatternPolicy/v1r3p3
GLASTPATTERNPOLICYVERSION=v1r3p3
GlastPatternPolicy_cmtpath=/nfs/farm/g/glast/u10/builds/GlastRelease/gr-v9r10-build
GlastPatternPolicy_offset=GlastPolicy
GlastPatternPolicy_project=Project1
GlastMain=${GLASTPOLICYROOT}/src/GlastMain.cxx
TestGlastMain=${GLASTPOLICYROOT}/src/TestGlastMain.cxx
GlastCppPolicy_tag=$(tag)
GLASTCPPPOLICYROOT=/nfs/farm/g/glast/u10/builds/GlastRelease/gr-v9r10-build/GlastPolicy/GlastCppPolicy/v1r6
GlastCppPolicy_root=/nfs/farm/g/glast/u10/builds/GlastRelease/gr-v9r10-build/GlastPolicy/GlastCppPolicy/v1r6
GLASTCPPPOLICYVERSION=v1r6
GlastCppPolicy_cmtpath=/nfs/farm/g/glast/u10/builds/GlastRelease/gr-v9r10-build
GlastCppPolicy_offset=GlastPolicy
GlastCppPolicy_project=Project1
BINDIR=rh9_gcc32
cppoptions=$(cppdebugflags_s)
cppdebugflags_s=-g
cppoptimized_s=-O2
cppprofiled_s=-pg
linkdebug=-g 
makeLinkMap=-Wl,-Map,Linux.map
componentshr_linkopts=-fPIC  -ldl 
libraryshr_linkopts=-fPIC -ldl 
geometry_linkopts=-L${geometry_cmtpath}/lib -lgeometry 
idents_tag=$(tag)
IDENTSROOT=/nfs/farm/g/glast/u10/builds/GlastRelease/gr-v9r10-build/idents/v2r18
idents_root=/nfs/farm/g/glast/u10/builds/GlastRelease/gr-v9r10-build/idents/v2r18
IDENTSVERSION=v2r18
idents_cmtpath=/nfs/farm/g/glast/u10/builds/GlastRelease/gr-v9r10-build
idents_project=Project1
facilities_tag=$(tag)
FACILITIESROOT=/nfs/farm/g/glast/u10/builds/GlastRelease/gr-v9r10-build/facilities/v2r12p4
facilities_root=/nfs/farm/g/glast/u10/builds/GlastRelease/gr-v9r10-build/facilities/v2r12p4
FACILITIESVERSION=v2r12p4
facilities_cmtpath=/nfs/farm/g/glast/u10/builds/GlastRelease/gr-v9r10-build
facilities_project=Project1
facilities_linkopts=-L${facilities_cmtpath}/lib -lfacilities 
idents_linkopts=-L${idents_cmtpath}/lib -lidents 
Event_linkopts=-L${Event_cmtpath}/lib -lEvent 
Event_shlibflags=$(libraryshr_linkopts)
Event_stamps=${EVENTROOT}/${BINDIR}/Event.stamp 
OnboardFilterTds_linkopts=-L${OnboardFilterTds_cmtpath}/lib -lOnboardFilterTds 
OnboardFilterTds_shlibflags=$(libraryshr_linkopts)
OnboardFilterTds_stamps=${ONBOARDFILTERTDSROOT}/${BINDIR}/OnboardFilterTds.stamp 
GlastPatternPolicyDir=${GLASTPATTERNPOLICYROOT}/${BINDIR}
GlastPolicyDir=${GLASTPOLICYROOT}/${BINDIR}
geometryDir=${GEOMETRYROOT}/${BINDIR}
facilitiesDir=${FACILITIESROOT}/${BINDIR}
identsDir=${IDENTSROOT}/${BINDIR}
EventDir=${EVENTROOT}/${BINDIR}
OnboardFilterTdsDir=${ONBOARDFILTERTDSROOT}/${BINDIR}
tag=rh9_gcc32
package=OnboardFilterTds
version=v0
PACKAGE_ROOT=$(ONBOARDFILTERTDSROOT)
srcdir=../src
bin=../$(OnboardFilterTds_tag)/
javabin=../classes/
mgrdir=cmt
project=Project1
use_requirements=requirements $(CMTROOT)/mgr/requirements $(OBFROOT)/cmt/requirements $(EVENTROOT)/cmt/requirements $(GAUDIINTERFACEROOT)/cmt/requirements $(EXTERNALLIBSROOT)/cmt/requirements $(ROOTROOT)/cmt/requirements $(GEOMETRYROOT)/mgr/requirements $(CLHEPROOT)/cmt/requirements $(CPPUNITROOT)/cmt/requirements $(XMLEXTROOT)/cmt/requirements $(IDENTSROOT)/cmt/requirements $(FACILITIESROOT)/cmt/requirements $(GLASTPOLICYROOT)/cmt/requirements $(GLASTPATTERNPOLICYROOT)/cmt/requirements $(GLASTCPPPOLICYROOT)/cmt/requirements 
use_includes= $(ppcmd)"$(obf_root)/src" $(ppcmd)"$(EVENTROOT)" $(ppcmd)"$(GaudiInterface_DIR)/include" $(ppcmd)"$(ROOT_PATH)/include" $(ppcmd)"$(GEOMETRYROOT)" $(ppcmd)"$(CLHEPBASE)/include" $(ppcmd)"${cppunit_DIR}/include" $(ppcmd)"$(XMLEXT_DIR)/include" $(ppcmd)"$(IDENTSROOT)" $(ppcmd)"$(FACILITIESROOT)" 
use_fincludes= $(use_includes)
use_stamps= $(OnboardFilterTds_stamps)  $(obf_stamps)  $(Event_stamps)  $(GaudiInterface_stamps)  $(ExternalLibs_stamps)  $(ROOT_stamps)  $(geometry_stamps)  $(CLHEP_stamps)  $(cppunit_stamps)  $(XMLEXT_stamps)  $(idents_stamps)  $(facilities_stamps)  $(GlastPolicy_stamps)  $(GlastPatternPolicy_stamps)  $(GlastCppPolicy_stamps) 
use_cflags=  $(OnboardFilterTds_cflags)  $(obf_cflags)  $(Event_cflags)  $(GaudiInterface_cflags)  $(ExternalLibs_cflags)  $(ROOT_cflags)  $(geometry_cflags)  $(CLHEP_cflags)  $(cppunit_cflags)  $(XMLEXT_cflags)  $(idents_cflags)  $(facilities_cflags)  $(GlastPolicy_cflags) 
use_pp_cflags=  $(OnboardFilterTds_pp_cflags)  $(obf_pp_cflags)  $(Event_pp_cflags)  $(GaudiInterface_pp_cflags)  $(ExternalLibs_pp_cflags)  $(ROOT_pp_cflags)  $(geometry_pp_cflags)  $(CLHEP_pp_cflags)  $(cppunit_pp_cflags)  $(XMLEXT_pp_cflags)  $(idents_pp_cflags)  $(facilities_pp_cflags)  $(GlastPolicy_pp_cflags) 
use_cppflags=  $(OnboardFilterTds_cppflags)  $(obf_cppflags)  $(Event_cppflags)  $(GaudiInterface_cppflags)  $(ExternalLibs_cppflags)  $(ROOT_cppflags)  $(geometry_cppflags)  $(CLHEP_cppflags)  $(cppunit_cppflags)  $(XMLEXT_cppflags)  $(idents_cppflags)  $(facilities_cppflags)  $(GlastPolicy_cppflags) 
use_pp_cppflags=  $(OnboardFilterTds_pp_cppflags)  $(obf_pp_cppflags)  $(Event_pp_cppflags)  $(GaudiInterface_pp_cppflags)  $(ExternalLibs_pp_cppflags)  $(ROOT_pp_cppflags)  $(geometry_pp_cppflags)  $(CLHEP_pp_cppflags)  $(cppunit_pp_cppflags)  $(XMLEXT_pp_cppflags)  $(idents_pp_cppflags)  $(facilities_pp_cppflags)  $(GlastPolicy_pp_cppflags) 
use_fflags=  $(OnboardFilterTds_fflags)  $(obf_fflags)  $(Event_fflags)  $(GaudiInterface_fflags)  $(ExternalLibs_fflags)  $(ROOT_fflags)  $(geometry_fflags)  $(CLHEP_fflags)  $(cppunit_fflags)  $(XMLEXT_fflags)  $(idents_fflags)  $(facilities_fflags)  $(GlastPolicy_fflags) 
use_pp_fflags=  $(OnboardFilterTds_pp_fflags)  $(obf_pp_fflags)  $(Event_pp_fflags)  $(GaudiInterface_pp_fflags)  $(ExternalLibs_pp_fflags)  $(ROOT_pp_fflags)  $(geometry_pp_fflags)  $(CLHEP_pp_fflags)  $(cppunit_pp_fflags)  $(XMLEXT_pp_fflags)  $(idents_pp_fflags)  $(facilities_pp_fflags)  $(GlastPolicy_pp_fflags) 
use_linkopts= $(cmt_installarea_linkopts)   $(OnboardFilterTds_linkopts)  $(obf_linkopts)  $(Event_linkopts)  $(GaudiInterface_linkopts)  $(ExternalLibs_linkopts)  $(ROOT_linkopts)  $(geometry_linkopts)  $(CLHEP_linkopts)  $(cppunit_linkopts)  $(XMLEXT_linkopts)  $(idents_linkopts)  $(facilities_linkopts)  $(GlastPolicy_linkopts) 
use_libraries= $(obf_libraries)  $(Event_libraries)  $(GaudiInterface_libraries)  $(ExternalLibs_libraries)  $(ROOT_libraries)  $(geometry_libraries)  $(CLHEP_libraries)  $(cppunit_libraries)  $(XMLEXT_libraries)  $(idents_libraries)  $(facilities_libraries)  $(GlastPolicy_libraries)  $(GlastPatternPolicy_libraries)  $(GlastCppPolicy_libraries) 
fincludes= $(includes)
OnboardFilterTds_use_linkopts=  $(OnboardFilterTds_linkopts)  $(obf_linkopts)  $(Event_linkopts)  $(GaudiInterface_linkopts)  $(ExternalLibs_linkopts)  $(ROOT_linkopts)  $(geometry_linkopts)  $(CLHEP_linkopts)  $(cppunit_linkopts)  $(XMLEXT_linkopts)  $(idents_linkopts)  $(facilities_linkopts)  $(GlastPolicy_linkopts) 
test_OnboardFilterTds_use_linkopts=  $(OnboardFilterTds_linkopts)  $(obf_linkopts)  $(Event_linkopts)  $(GaudiInterface_linkopts)  $(ExternalLibs_linkopts)  $(ROOT_linkopts)  $(geometry_linkopts)  $(CLHEP_linkopts)  $(cppunit_linkopts)  $(XMLEXT_linkopts)  $(idents_linkopts)  $(facilities_linkopts)  $(GlastPolicy_linkopts) 
constituents= OnboardFilterTds test_OnboardFilterTds 
all_constituents= $(constituents)
constituentsclean= test_OnboardFilterTdsclean OnboardFilterTdsclean 
all_constituentsclean= $(constituentsclean)
cmt_installarea_paths=$(cmt_installarea_prefix)/$(tag)/bin $(cmt_installarea_prefix)/$(tag)/lib $(cmt_installarea_prefix)/share/lib $(cmt_installarea_prefix)/share/bin
