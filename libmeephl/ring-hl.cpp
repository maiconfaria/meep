/***************************************************************/
/* C++ port of meep/examples/ring.ctl, using the               */
/* "high-level" meep C++ interface stack, which consists of    */
/* libmeep_hl + libmeep_geom + libctlgeom + libmeep            */
/***************************************************************/
/*
; Calculating 2d ring-resonator modes, from the Meep tutorial.
*/

#include <stdio.h>
#include <stdlib.h>
#include <complex>

#include "meephl.hpp"

using namespace meep;
using namespace meephl;

typedef std::complex<double> cdouble;

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[])
{
  meep_session S(argc, argv);

  double n=3.4;     // index of waveguide
  double w=1.0;     // width of waveguide
  double r=1.0;     // inner radius of ring

  double pad=4;     // padding between waveguide and edge of PML
  double dpml=2;    // thickness of PML

  double sxy        = 2.0*(r+w+pad+dpml);  // cell size

  // (set-param! resolution 10)
  // (set! geometry-lattice (make lattice (size sxy sxy no-size)))
  // (set! symmetries (list (make mirror-sym (direction Y))))
  S.set_resolution(10.0);
  S.set_geometry_lattice(sxy, sxy);
  S.add_mirror_symmetry(Y);

  // ; Create a ring waveguide by two overlapping cylinders - later objects
  // ; take precedence over earlier objects, so we put the outer cylinder first.
  // ; and the inner (air) cylinder second.
  // (set! geometry (list
  //	(make cylinder (center 0 0) (height infinity)
  //		(radius (+ r w)) (material (make dielectric (index n))))
  // 	(make cylinder (center 0 0) (height infinity)
  // 		(radius r) (material air))))
  material_type dielectric = meep_geom::make_dielectric(n*n);
  material_type vacuum     = meep_geom::vacuum;
  S.add_object(dielectric, "CYLINDER --height %e --radius %e", HUGE_VAL, r);
  S.add_object(vacuum,     "CYLINDER --height %e --radius %e", HUGE_VAL, r+w);

  // (set! sources (list
  //              (make source
  //                (src (make gaussian-src (frequency fcen) (fwidth df)))
  //                (component Ez) (center (+ r 0.1) 0))))
  double fcen = 0.15;  // ; pulse center frequency
  double df   = 0.1;   // ; df
  S.add_gaussian_src(fcen, df, Ez, vec(r+0.1, 0.0), vec(0.0, 0.0));

  // (run-sources+ 300 
  // 	(at-beginning output-epsilon)
  // 	(after-sources (harminv Ez (vector3 (+ r 0.1)) fcen df)))
  S.add_output(Dielectric, AT_BEGINNING);
  S.run_sources_plus(300.0);

  // (run-until (/ 1 fcen) (at-every (/ 1 fcen 20) output-efield-z))
  S.add_output(Ez, AT_EVERY, fcen/20.0);
  S.run_until(1.0/fcen);

  // success if we made it here
  exit(0);

}
