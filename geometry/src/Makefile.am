## Process this file with automake to produce Makefile.in
#
# $Id$
#
# Author: Paul_Kunz@slac.stanford.edu
#

# The following set,  else would follow GNU conventions
AUTOMAKE_OPTIONS = foreign

noinst_LTLIBRARIES = libgeometry.la

## leave `Hype.cxx' out since many functions are not implemented
libgeometry_la_SOURCES =	\
	Box.cxx			\
	Cone.cxx		\
	Conic.cxx		\
	Cons.cxx		\
	CoordTransform.cxx	\
	Cylinder.cxx		\
	GeomObject.cxx		\
	Helix.cxx		\
	Hexahedron.cxx		\
	Hyperbolic.cxx		\
	Intersection.cxx	\
	Plane.cxx		\
	Point.cxx		\
	Ray.cxx			\
	Shape.cxx		\
	Sphe.cxx		\
	Sphere.cxx		\
	Surface.cxx		\
	Track.cxx		\
	Trap.cxx		\
	Trd1.cxx		\
	Trd2.cxx		\
	Tube.cxx		\
	Tubs.cxx		\
	Vector.cxx		\
	Volume.cxx		\
	Wedge.cxx

INCLUDES = -I$(top_srcdir)		\
	-I$(top_srcdir)/../CLHEP	\
	-I$(glast_prefix)/include
