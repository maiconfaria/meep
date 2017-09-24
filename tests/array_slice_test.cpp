#include <math.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <meep.hpp>
#include "meep_internals.hpp"
#include "config.h"

using namespace meep;
using namespace std;

typedef complex<double> cdouble;

const double xsize = 2.0;
const double ysize = 2.0;
const double zsize = 0.6;

const double r = 0.5;
const double eps_k = 2*pi / 1.0;

double funky_eps_2d(const vec &p_) {
  vec p = p_ - vec(xsize / 2, ysize / 2);
  if (fabs(p & p) < r * r)
    return 1.0;
  return 2.0 + cos(p.x() * eps_k) * cos(p.y() * eps_k);
}

double funky_eps_3d(const vec &p_) {
  vec p = p_ - vec(xsize / 2, ysize / 2, zsize / 2);
  if (fabs(p & p) < r * r)
    return 1.0;
  return 2.0 + cos(p.x() * eps_k) * cos(p.y() * eps_k) * cos(p.z() * eps_k);
}

symmetry make_identity(const grid_volume &gv)
{
  (void) gv; // unused
  return identity();
}

symmetry make_mirrorx(const grid_volume &gv)
{
  return mirror(X, gv);
}

symmetry make_mirrory(const grid_volume &gv)
{
  return mirror(Y, gv);
}

symmetry make_mirrorxy(const grid_volume &gv)
{
  return mirror(X, gv) + mirror(Y, gv);
}

symmetry make_rotate4z(const grid_volume &gv)
{
  return rotate4(Z, gv);
}

typedef symmetry (*symfunc)(const grid_volume &);

/***************************************************************/
/* compare two real-valued numbers *****************************/
/***************************************************************/
const double tol = sizeof(realnum) == sizeof(float) ? 1e-4 : 1e-8;

double compare(double a, double b, const char *nam, int i0,int i1,int i2) {
  if (fabs(a-b) > tol*tol + fabs(b) * tol || b != b) {
    master_printf("%g vs. %g differs by\t%g\n", a, b, fabs(a-b));
    master_printf("This gives a fractional error of %g\n", fabs(a-b)/fabs(b));
    abort("Error in %s at (%d,%d,%d)\n", nam, i0,i1,i2);
  }
  return fabs(a-b);
}

/***************************************************************/
/* compare two complex-valued numbers **************************/
/***************************************************************/
double compare(cdouble a, cdouble b, const char *nam, int i0,int i1,int i2) {
  if (abs(a-b) > tol*tol + abs(b) * tol || b != b) {
    master_printf("{%g,%g} vs. {%g,%g} differs by\t%g\n", real(a), imag(a), real(b), imag(b), abs(a-b));
    master_printf("This gives a fractional error of %g\n", abs(a-b)/abs(b));
    abort("Error in %s at (%d,%d,%d)\n", nam, i0,i1,i2);
  }
  return abs(a-b);
}

double get_reim(complex<double> x, int reim)
{
  return reim ? imag(x) : real(x);
}

/***************************************************************/
/* modeled after check_2d in h5test.cpp  ***********************/
/***************************************************************/
bool check_2d(double eps(const vec &), double a, int splitting, symfunc Sf,
	      double kx, double ky, component src_c, int eval_c,
	      volume where, bool real_fields, int expected_rank,
	      const char *casename) {
  
  (void) expected_rank;

  /***************************************************************/
  /* initialize structure, fields, sources and run calculation   */
  /***************************************************************/
  const grid_volume gv = vol2d(xsize, ysize, a);
  structure s(gv, eps, no_pml(), Sf(gv), splitting);
  fields f(&s);

  f.use_bloch(X, real_fields ? 0.0 : kx);
  f.use_bloch(Y, real_fields ? 0.0 : ky);

  if (real_fields) f.use_real_fields();
  f.add_point_source(src_c, 0.3, 2.0, 0.0, 1.0, gv.center(), 1.0, 1);

  if (eval_c >= int(Dielectric)) real_fields = true;

  while (f.time() <= 3.0 && !interrupt)
    f.step();

  bool has_imag = !(f.is_real) && (eval_c != Dielectric) && (eval_c != Permeability);
printf("has_imag=%s\n",has_imag ? "true" : "false");

  /***************************************************************/
  /* fetch array slice using get_array_slice *********************/
  /***************************************************************/
  int dims[3];
  int rank=f.get_array_slice_dimensions(where, dims);
  if (rank!=2)
   abort("case %s: rank=%i (should be %i)",casename,rank,2);
  double *slice; 
  cdouble *zslice;
  component c=component(eval_c);
  if (has_imag)
   zslice=f.get_complex_array_slice(where, c);
  else
   slice=f.get_array_slice(where, c);

  /***************************************************************/
  /* test each element of the array slice against the            */
  /* corresponding element in the fields array                   */
  /***************************************************************/

  // compute corner coordinate of slice
  vec loc0(where.get_min_corner());
  ivec iloc0(gv.dim);
  LOOP_OVER_DIRECTIONS(gv.dim, d) {
    iloc0.set_direction(d, 1+2*int(floor(loc0.in_direction(d)*a-.5)));
    if (where.in_direction(d) == 0.0 &&
	1. - where.in_direction_min(d)*a + 0.5*iloc0.in_direction(d)
	<= 1. + where.in_direction_max(d)*a - 0.5*(iloc0.in_direction(d)+2))
      iloc0.set_direction(d, iloc0.in_direction(d) + 2); // snap to grid
  }
  loc0 = gv[iloc0];

  double data_min = meep::infinity, data_max = -meep::infinity;
  double err_max = 0;
  for (int reim = 0; reim < (real_fields ? 1 : 2); ++reim) {
    vec loc(loc0.dim);
    for (int i0 = 0; i0 < dims[0]; ++i0) {
      for (int i1 = 0; i1 < dims[1]; ++i1) {
	loc.set_direction(X, loc0.in_direction(X) + i0 * gv.inva);
	loc.set_direction(Y, loc0.in_direction(Y) + i1 * gv.inva);
	int idx = i0 * dims[1] + i1;

	/* Ugh, for rotational symmetries (which mix up components etc.),
	   we can't guarantee that a component is *exactly* the
	   same as its rotated version, and we don't know which one
	   was written to the file (HR 20170921 file->slice). */
	int cs = eval_c;
	complex<double> ph = 1.0;
	double correct_value=get_reim(f.get_field(eval_c,loc),reim);
	double   slice_value=has_imag ? get_reim(zslice[idx],reim)
                                      : slice[idx];
	double diff = fabs(correct_value-slice_value);
	for (int sn = 1; sn < f.S.multiplicity(); ++sn) {
	  vec loc2(f.S.transform(loc, sn));
	  int cs2 = f.S.transform(eval_c, sn);
	  complex<double> ph2 = f.S.phase_shift(cs2, -sn);
	  double correct_value2=get_reim(f.get_field(cs2,loc2)*ph2,reim);
	  double diff2 = fabs(correct_value2-slice_value);

	  if (diff2 < diff) {
	    loc = loc2;
	    cs = cs2;
	    ph = ph2;
	    diff = diff2;
	  }
	};

	double err = compare(slice_value,
			     get_reim(f.get_field(cs, loc) * ph, reim),
			     casename,i0,i1,0);
	err_max = max(err, err_max);
	data_min = min(data_min, slice_value);
	data_max = max(data_max, slice_value);
      }
    }
  };

  if (has_imag)
   delete zslice;
  else
   delete slice;

  master_printf("Passed case %s, cmpt %s (%g..%g), err=%g\n", casename, component_name(c),
		data_min, data_max,
		err_max / max(fabs(data_min), fabs(data_max)));

  return true;
}

#if 0
bool check_3d(double eps(const vec &), double a, int splitting, symfunc Sf,
	      component src_c, int eval_c,
	      volume file_gv,
	      bool real_fields, int expected_rank,
	      const char *name) {
  const grid_volume gv = vol3d(xsize, ysize, zsize, a);
  structure s(gv, eps, no_pml(), Sf(gv), splitting);
  fields f(&s);

  if (real_fields) f.use_real_fields();
  f.add_point_source(src_c, 0.3, 2.0, 0.0, 1.0, gv.center(), 1.0, 1);

  if (eval_c >= Dielectric) real_fields = true;

  while (f.time() <= 3.0 && !interrupt)
    f.step();

  h5file *file = f.open_h5file(name);
  if (is_derived(eval_c))
    f.output_hdf5(derived_component(eval_c), file_gv, file);
  else
    f.output_hdf5(component(eval_c), file_gv, file);

  file->write("stringtest", "Hello, world!\n");

  delete file;
  all_wait();
  sync();
  file = f.open_h5file(name, h5file::READONLY);

  char *str = file->read("stringtest");
  if (strcmp(str, "Hello, world!\n"))
       abort("Failed to read back string test from %s...", name);

  // compute corner coordinate of file data
  vec loc0(file_gv.get_min_corner());
  ivec iloc0(gv.dim);
  LOOP_OVER_DIRECTIONS(gv.dim, d) {
    iloc0.set_direction(d, 1+2*int(floor(loc0.in_direction(d)*a-.5)));
    if (file_gv.in_direction(d) == 0.0 &&
	1. - file_gv.in_direction_min(d)*a + 0.5*iloc0.in_direction(d)
	<= 1. + file_gv.in_direction_max(d)*a - 0.5*(iloc0.in_direction(d)+2))
      iloc0.set_direction(d, iloc0.in_direction(d) + 2); // snap to grid
  }
  loc0 = gv[iloc0];

  double data_min = meep::infinity, data_max = -meep::infinity;
  double err_max = 0;
  for (int reim = 0; reim < (real_fields ? 1 : 2); ++reim) {
    int rank, dims[3] = {1, 1, 1};

    char dataname[256];
    snprintf(dataname, 256, "%s%s", component_name(eval_c),
	     reim ? ".i" : (real_fields ? "" : ".r"));

    realnum *h5data = file->read(dataname, &rank, dims, 3);
    file->prevent_deadlock(); // hackery
    if (!h5data)
	 abort("failed to read dataset %s:%s\n", name, dataname);
    if (rank != expected_rank)
	 abort("incorrect rank (%d instead of %d) in %s:%s\n",
	       rank, expected_rank, name, dataname);
    vec loc(loc0.dim);
    for (int i0 = 0; i0 < dims[0]; ++i0) {
      for (int i1 = 0; i1 < dims[1]; ++i1) {
	for (int i2 = 0; i2 < dims[2]; ++i2) {
	  loc.set_direction(X, loc0.in_direction(X) + i0 * gv.inva);
	  loc.set_direction(Y, loc0.in_direction(Y) + i1 * gv.inva);
	  loc.set_direction(Z, loc0.in_direction(Z) + i2 * gv.inva);
	  int idx = (i0 * dims[1] + i1) * dims[2] + i2;
	  
	  /* Ugh, for rotational symmetries (which mix up components etc.),
	     we can't guarantee that a component is *exactly* the
	     same as its rotated version, and we don't know which one
	     was written to the file. */
	  int cs = eval_c;
	  complex<double> ph = 1.0;
	  double diff = fabs(get_reim(f.get_field(eval_c, loc), reim) -
			     h5data[idx]);
	  for (int sn = 1; sn < f.S.multiplicity(); ++sn) {
	    vec loc2(f.S.transform(loc, sn));
	    int cs2 = f.S.transform(eval_c, sn);
	    complex<double> ph2 = f.S.phase_shift(cs2, -sn);
	    double diff2 = fabs(get_reim(f.get_field(cs2, loc2)*ph2, reim) -
				h5data[idx]);
	    if (diff2 < diff) {
	      loc = loc2;
	      cs = cs2;
	      ph = ph2;
	      diff = diff2;
	    }
	  }
	  
	  double err = compare(h5data[idx],
			       get_reim(f.get_field(cs, loc)*ph,reim),
			       name, i0,i1,i2);
	  err_max = max(err, err_max);
	  data_min = min(data_min, h5data[idx]);
	  data_max = max(data_max, h5data[idx]);
	}
      }
    }
    delete[] h5data;
  }

  //file->remove();
  delete file;

  master_printf("Passed %s (%g..%g), err=%g\n", name,
		data_min, data_max,
		err_max / (max(fabs(data_min), fabs(data_max)) + 1e-16));

  return 1;
}
#endif

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char **argv)
{
  initialize mpi(argc, argv);

  /***************************************************************/
  /* allow user to say e.g. --single2dtest 0 1 0 1 1             */
  /* to run just a single one of the 2D tests below              */
  /***************************************************************/
  int which2dtest[5];
  bool single2dtest=false;
  for(int narg=1; narg<argc; narg++)
   if ( !strcmp(argv[narg],"--single2dtest") )
    { single2dtest=true;
      if (argc<narg+6) abort("invalid command-line syntax");
      sscanf(argv[narg+1],"%i",which2dtest+0);
      sscanf(argv[narg+2],"%i",which2dtest+1);
      sscanf(argv[narg+3],"%i",which2dtest+2);
      sscanf(argv[narg+4],"%i",which2dtest+3);
      sscanf(argv[narg+5],"%i",which2dtest+4);
      printf("Running single 2d test (%i,%i,%i,%i,%i)\n",
              which2dtest[0], which2dtest[1], which2dtest[2],
              which2dtest[3], which2dtest[4]);
    };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  const double a = 10.0;
  int chances;
  quiet = true;
  const double pad1 = 0.314159, pad2 = 0.27183, pad3 = 0.14142;

  volume gv_2d[4] = {
       volume(vec(pad1,pad2), vec(xsize-pad2,ysize-pad1)),
       volume(vec(-pad1,-pad2), vec(2*xsize-pad2,2*ysize-pad1)),
       volume(vec(pad1,pad2), vec(xsize-pad2,pad2)),
       volume(vec(pad1,pad2), vec(pad1,pad2)),
  };
  char gv_2d_name[4][20] = {"plane", "plane-supercell", "line", "point"};
  int gv_2d_rank[4] = {2,2,1,0};
  //int tm_c[5] = {Dielectric, Ez, Hy, Sx, D_EnergyDensity};
  int tm_c[5] = {Ez, Dielectric, Hy, Sx, D_EnergyDensity};
  symfunc Sf2[5] = {make_identity, make_mirrorx, make_mirrory, make_mirrorxy,
		   make_rotate4z};
  char Sf2_name[5][32] = {"identity", "mirrorx", "mirrory", "mirrorxy",
			 "rotate4z"};
  double Sf2_kx[5] = {0.3, 0, 0.3, 0, 0};
  double Sf2_ky[5] = {0.2, 0.2, 0, 0, 0};

  /* 
     if the user didn't specify --single2dtest, then... 
     this test takes too long, so only do 1/chances of the cases,
     randomly selected 
  */
  //srand(314159); /* deterministic "rand" */
  srand(time(0)); // HR 20170921 make it more interesting
  chances = argc > 1 ? atoi(argv[1]) : 5;
  for (int iS = 0; iS < 5; ++iS)
    for (int splitting = 0; splitting < 5; ++splitting)
      for (int igv = 0; igv < 4; ++igv)
	for (int ic = 0; ic < 5; ++ic)
	  for (int use_real = 1; use_real >= 0; --use_real)
           { 
             bool skip;
	     if ( single2dtest )
	      skip = (     which2dtest[0]!=iS
	               ||  which2dtest[1]!=splitting
	               ||  which2dtest[2]!=igv
	               ||  which2dtest[3]!=ic
	               ||  which2dtest[4]!=use_real
	             );
	     else
	      skip = (broadcast(0, rand()) % chances != 0);

	     if (skip) continue;

	     char casename[1024];
	     snprintf(casename, 1024, "check_2d_tm_%s_%d_%s_%s%s",
		      Sf2_name[iS], splitting, gv_2d_name[igv],
		      component_name(tm_c[ic]), use_real ? "_r" : "");
	     master_printf("Checking case %s...\n", casename);
	     if (!check_2d(funky_eps_2d, a, splitting,
			   Sf2[iS], Sf2_kx[iS], Sf2_ky[iS],
			   Ez, tm_c[ic], gv_2d[igv],
			   use_real, gv_2d_rank[igv], casename)
                ) return 1;
	   };
#if 0
  for (int iS = 0; iS < 3; ++iS)
    for (int splitting = 0; splitting < 5; splitting += 3)
      for (int igv = 0; igv < 4; ++igv) {
	for (int ic = 0; ic < 1; ++ic)
	  if (broadcast(0, rand()) % chances == 0) {
	    bool use_real = true;
	    char name[1024];
	    snprintf(name, 1024, "check_3d_ezsrc_%s_%d_%s_%s%s", Sf3_name[iS],
		     splitting, gv_3d_name[igv], component_name(c3d[ic]),
		     use_real ? "_r" : "");
	    master_printf("Checking %s...\n", name);
	    if (!check_3d(funky_eps_3d, a, splitting, Sf3[iS], Ez, c3d[ic],
	    		  gv_3d[igv], use_real, gv_3d_rank[igv], name))
	      return 1;
	  }
      }
#endif
  return 0;
}