LIBCTLGEOM   = -lctlgeom
LIBMEEP      = $(top_builddir)/src/libmeep.la
LIBMEEPGEOM  = $(top_builddir)/libmeepgeom/libmeepgeom.la
MEEPLIBS     = $(LIBMEEPGEOM) $(LIBCTLGEOM) $(LIBMEEP)

LIBSCUFF     = -lscuff
SCUFFLIBS    = $(LIBSCUFF)

lib_LTLIBRARIES = libmeephl.la
pkginclude_HEADERS = meephl.hpp

AM_CPPFLAGS = -I$(top_srcdir)/src                   \
              -I$(top_srcdir)/libmeepgeom           \
              -I$(SCUFF)/include/scuff-em

libmeephl_la_SOURCES = \
 meephl.cpp meephl.hpp

libmeephl_la_LDFLAGS = -version-info @SHARED_VERSION_INFO@ -L$(SCUFF)/lib

check_PROGRAMS  = ring-hl

ring_hl_SOURCES = ring-hl.cpp
ring_hl_LDADD   = libmeephl.la $(MEEPLIBS) $(SCUFFLIBS)

noinst_PROGRAMS = ring-hl

TESTS           = ring-hl

#SWIG_SRC = meephl.i
#EXTRA_DIST = $(SWIG_SRC) 
#BUILT_SOURCES = meephl.cpp meephl.py
#
#_meephl_la_SOURCES = meephl.cpp
#_meephl_la_LIBADD = $(MEEPLIBS) $(SCUFFLIBS) $(LIBS)
#_meephl_la_LDFLAGS = -module -version-info @SHARED_VERSION_INFO@
#_meephl_la_CPPFLAGS = $(PYTHON_INCLUDES) $(AM_CPPFLAGS)
#
#python_PYTHON = meephl.py
#pyexec_LTLIBRARIES = _meephl.la
#
#if MAINTAINER_MODE
#
#meephl.cpp meephl.py: $(SWIG_SRC) meephl.hpp
#	swig -Wall $(AM_CPPFLAGS) -outdir $(builddir) -c++ -python -o $(builddir)/meephl.cpp $(srcdir)/meephl.i
#
#endif
