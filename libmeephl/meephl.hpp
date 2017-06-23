/* Copyright (C) 2005-2015 Massachusetts Institute of Technology
%
%  This program is free software; you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation; either version 2, or (at your option)
%  any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program; if not, write to the Free Software Foundation,
%  Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
*/
#ifndef MEEPHL_H
#define MEEPHL_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>

#include <meep.hpp>
#include <meep/vec.hpp>
#include <ctlgeom.h>
#include <meep/meepgeom.hpp>

#define NO_SIZE 1.0e-20

namespace meephl {

// prototypes for user-supplied routines

// permittivity function
typedef double (*eps_func)(const meep::vec &x, void *user_data);

// step function
typedef void (*step_func)(meep::fields *f, void *user_data);

// condition function (return true if timestepping should end)
typedef bool (*cond_func)(meep::fields *f, void *user_data);

/***************************************************************/
/***************************************************************/
/***************************************************************/
class user_material_function : public meep::material_function
{
  // class methods 
public:
  user_material_function(eps_func func, void *user_data);
  virtual double eps(const meep::vec &r); 

  // class fields
//private:
  eps_func func;
  void *user_data;
};

// when to run step functions or output components
typedef enum 
 { AT_BEGINNING, AT_END, AT_TIME, AT_EVERY } runtime;

/***************************************************************/
/* data structures used internally within meep_session **********/
/***************************************************************/
typedef struct step_func_data
 {
   step_func func;
   void *user_data;
   runtime when;
   double T;
   double TNext;
   step_func_data *next;

 } step_func_data;

/***************************************************************/
/***************************************************************/
/***************************************************************/
typedef struct source_data  
 {
   double fcen;
   double df;
   meep::component c;
   meep::vec center;
   meep::vec size;
   source_data *next;

 } source_data; 
/***************************************************************/
/***************************************************************/
/***************************************************************/
class meep_session
 {
    public:
       meep_session(int argc, char *argv[]);
      ~meep_session();

  /***************************************************************/
  /*- 1) routines to set up the geometry before timestepping     */
  /***************************************************************/

       /*-----------------------------------*/
       /*-- configure geometric parameters  */
       /*-----------------------------------*/
       void set_resolution(double res);
       void set_geometry_lattice(int nz);
       void set_geometry_lattice(int nx, int ny);
       void set_geometry_lattice(int nx, int ny, int nz);
       void set_geometry_cylindrical(int nr, int nz=NO_SIZE);

       void add_mirror_symmetry(meep::direction d);
       void add_rotate2_symmetry(meep::direction d);
       void add_rotate4_symmetry(meep::direction d);

       void add_pml(double thickness);
       void add_pml(double thickness, meep::direction d);
       void add_pml(double thickness, meep::direction d, meep::boundary_side side);

       /*-----------------------------------------*/
       /*- define the material geometry, either   */
       /*- by providing a position-dependent      */
       /*- permittivity function or by specifying */
       /*- a list of geometric objects            */
       /*-----------------------------------------*/
       void set_eps_func(eps_func func, void *user_data);

       geometric_object make_object(const char *description, ...);

       geometric_object *add_object(material_type material,
                                    const char *description, ...);

       void add_object(material_type material,
                       geometric_object *o);

       geometric_object *transform_object(geometric_object *o,
                                          const char *transform_string,
                                          ...);

  /***************************************************************/
  /*- 2) routines to add sources                                 */
  /***************************************************************/
       void add_gaussian_src(double fcen, double df,
                             meep::component c,
                             meep::vec center, meep::vec size);

  /***************************************************************/
  /*- 3) routines to add step functions, either user-defined or  */
  /*-    built-in.                                               */
  /*-    in each case the former prototype adds to the internal  */
  /*-    default list of step functions.                         */
  /***************************************************************/
       step_func_data new_step_func_list();

       void add_step_func(step_func func, void *user_data,
                          runtime when, double T=0);

       void add_step_func(step_func_data **sf_list,
                          step_func func, void *user_data,
                          runtime when, double T=0);

       void add_output(meep::component c, 
                       runtime when, double T=0.0);

       void add_output(step_func_data **sf_list, meep::component c,
                       runtime when, double T=0.0);

       void clear_step_funcs();

  /***************************************************************/
  /** 4) timestepping routines.                                  */
  /***************************************************************/
       void run_until(step_func_data *step_func_list, cond_func func, void *user_data);

       void run_until(cond_func func, void *user_data);

       void run_until(step_func_data *step_func_list, double T);
       void run_until(double T);

       void run_sources(step_func_data *step_func_list);
       void run_sources();

       void run_sources_plus(step_func_data *step_func_list, double T);
       void run_sources_plus(double T);

       void run_until_fields_decayed(step_func_data *step_func_list,
                                     double dT, meep::component c,
                                     meep::vec pt, double decay_by);

       void run_until_fields_decayed(double dT, meep::component c,
                                     meep::vec pt, double decay_by);

  /*--------------------------------------------------------------*/
  /*- the remaining routines are intended for internal use and    */
  /*- would be marked 'private' if i cared about that distinction */
  /*--------------------------------------------------------------*/
//private:
       void init_structure();
       void init_fields();
       void call_step_funcs(step_func_data *sf_list, double TNow);
  
  /*--------------------------------------------------------------*/
  /*-  class data fields -----------------------------------------*/
  /*--------------------------------------------------------------*/
//private:
       double resolution;
       meep::grid_volume gv;
       meep::symmetry sym;
       meep::boundary_region br;
       user_material_function *the_material_function;
       geometric_object_list the_object_list;

       source_data *sources_to_add;

       step_func_data *step_func_list;

       meep::structure *the_structure;
       meep::fields    *the_fields;

  // static class variables
       static bool mpiInitialized;
 };

}; // namespace meephl

#endif // #ifndef MEEPHL_H
