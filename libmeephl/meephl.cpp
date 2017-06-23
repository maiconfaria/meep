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
%  GNU General Public License for more details.  %
%  You should have received a copy of the GNU General Public License
%  along with this program; if not, write to the Free Software Foundation,
%  Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
*/

#include "meephl.hpp"
#include <libhrutil.h>

namespace meephl {

double dummy_eps_func(const meep::vec &x, void *user_data)
{ return 1.0; }

/***************************************************************/
/* implementation of user_material_function class              */
/***************************************************************/
user_material_function::user_material_function(eps_func the_func, void *the_user_data)
{ func = the_func;
  user_data = the_user_data;
}

double user_material_function::eps(const meep::vec &r)
 { return func(r, user_data); }

/***************************************************************/
/* static class variables **************************************/
/***************************************************************/
bool meep_session::mpiInitialized=false;

/***************************************************************/
/***************************************************************/
/***************************************************************/
meep_session::meep_session(int argc, char *argv[])
{
  if (!mpiInitialized)
   { mpiInitialized=true;
     meep::initialize mpi(argc, argv);
   };
  
  //gv.center_origin();
  // default values for class variables
  resolution = 10.0;

  the_material_function = new user_material_function(dummy_eps_func, 0);
  the_object_list.num_items=0;
  the_object_list.items=0;
  
  sources_to_add=0;

  step_func_list=0;

  the_structure=0;
  the_fields=0;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
meep_session::~meep_session()
{
}

/*--------------------------------------------------*/
/*   routines that configure geometric parameters   */
/*--------------------------------------------------*/
void meep_session::set_resolution(double res)
{ resolution=res; }

void meep_session::set_geometry_lattice(int nz)
{ gv=meep::vol1d(nz, resolution); gv.center_origin(); }

void meep_session::set_geometry_lattice(int nx, int ny)
{ gv=meep::vol2d(nx, ny, resolution); gv.center_origin(); }

void meep_session::set_geometry_lattice(int nx, int ny, int nz)
{ gv=meep::vol3d(nx, ny, nz, resolution); gv.center_origin(); }

void meep_session::set_geometry_cylindrical(int nr, int nz)
{ gv=meep::volcyl(nr, nz, resolution); gv.center_origin(); }

/*--------------------------------------------------*/
/*- symmetries -------------------------------------*/
/*--------------------------------------------------*/
void meep_session::add_mirror_symmetry(meep::direction d)
{ sym = sym + meep::mirror(d, gv); }

void meep_session::add_rotate2_symmetry(meep::direction d)
{ sym = sym + meep::rotate2(d, gv); }

void meep_session::add_rotate4_symmetry(meep::direction d)
{ sym = sym + meep::rotate4(d, gv); }


/*--------------------------------------------------*/
/*- PML --------------------------------------------*/
/*--------------------------------------------------*/
void meep_session::add_pml(double thickness)
{ br = meep::pml(thickness); }

void meep_session::add_pml(double thickness, meep::direction d)
{ br = br + meep::pml(thickness, d); }

void meep_session::add_pml(double thickness, meep::direction d, meep::boundary_side side)
{ br = br + meep::pml(thickness, d, side); }

/*--------------------------------------------------*/
/*- position-dependent material specification,      */
/*- defined either by                               */
/*-  (a) user-supplied function                     */
/*- or                                              */
/*-  (b) user-supplied list of geometric_objects    */
/*--------------------------------------------------*/
void meep_session::set_eps_func(eps_func func, void *user_data)
{
  the_material_function=new user_material_function(func, user_data);
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void meep_session::add_gaussian_src(double fcen,
                                   double df,
                                   meep::component c,
                                   meep::vec center, 
                                   meep::vec size)
{
  source_data *sd = (source_data *)malloc(sizeof(*sd));
  sd->fcen   = fcen;
  sd->df     = df;
  sd->c      = c;
  sd->center = center;
  sd->size   = size;
  sd->next   = 0;

  if (sources_to_add==0)
   sources_to_add=sd->next;
  else
   { source_data *list_tail=sources_to_add;
     while(list_tail->next != 0) 
      list_tail=list_tail->next;
     list_tail->next = sd;
   };
}

/*--------------------------------------------------*/
/*--------------------------------------------------*/
/*--------------------------------------------------*/
void meep_session::init_structure()
{ 
  the_structure = new meep::structure(gv, *the_material_function, br, sym);

  if (the_object_list.num_items>0)
   meep_geom::set_materials_from_geometry(the_structure,
                                          the_object_list);
}

/*--------------------------------------------------*/
/*--------------------------------------------------*/
/*--------------------------------------------------*/
void meep_session::init_fields()
{
  if (!the_structure)
   init_structure();

  the_fields = new meep::fields(the_structure);
 
  source_data *fd=sources_to_add;
  while(fd)
   { meep::gaussian_src_time src(fd->fcen, fd->df);
     meep::volume v(fd->center, fd->size);
     the_fields->add_volume_source(fd->c, src, v);
     source_data *thisfd = fd;
     fd=fd->next;
     free(thisfd);
   };
  sources_to_add=0;

}

/*--------------------------------------------------*/
/*- step functions ---------------------------------*/
/*--------------------------------------------------*/
bool check_time(step_func_data *sfd, double TNow)
{ 
  double T=sfd->T, TNext=sfd->TNext;
  switch(sfd->when)
   { 
     case AT_BEGINNING:
       return TNow==-HUGE_VAL;

     case AT_TIME:
       return TNow>=T;

     case AT_EVERY:
       if ( TNow < TNext )
        return false;
       sfd->TNext += T;
       return true;

     case AT_END:
       return TNow==+HUGE_VAL;
   }; 

  return false;
}

void meep_session::call_step_funcs(step_func_data *sf_list, 
                                  double TNow)
{
  for(step_func_data *sfd=sf_list; sfd; sfd=sfd->next)
   if( check_time(sfd, TNow) )
    sfd->func(the_fields, sfd->user_data);
}

void meep_session::add_step_func(step_func_data **sf_list,
                                step_func func, void *user_data,
                                runtime when, double T)
{ 
  step_func_data *data = (step_func_data *)malloc(sizeof(*data));
  data->func      = func;
  data->user_data = user_data;
  data->when      = when;
  data->T         = T;
  data->TNext     = T;
  data->next      = 0;

  if ( *(sf_list)==0 )
   *sf_list = data;
  else
   { step_func_data *sfd = (*sf_list);
     while( sfd->next !=0 ) 
      sfd=sfd->next;
     sfd->next=data;
   };
}

void meep_session::add_step_func(step_func func, void *user_data,
                                runtime when, double T)
{ add_step_func(&step_func_list, func, user_data, when, T); }

void meep_session::clear_step_funcs()
{ step_func_list=0; /* TODO: properly deallocate nodes*/ }

/*--------------------------------------------------*/
/*- outputs (special case of step function)        -*/
/*--------------------------------------------------*/
typedef struct output_sf_data
 {
   meep::component c;

 } output_sf_data;

void output_sf(meep::fields *f, void *user_data)
{
  output_sf_data *data=(output_sf_data *)user_data;

  f->output_hdf5(data->c, f->total_volume());
}

void meep_session::add_output(step_func_data **sf_list,
                             meep::component c, runtime when,
                             double T)
{
  output_sf_data *data=(output_sf_data *)malloc(sizeof(*data));
  data->c=c;
  add_step_func(sf_list, output_sf, (void *)data, when, T);
}

void meep_session::add_output(meep::component c, runtime when,
                             double T)
{ add_output(&step_func_list, c, when, T); }

/*--------------------------------------------------*/
/* timestepping ------------------------------------*/
/*--------------------------------------------------*/
void meep_session::run_until(step_func_data *sf_list,
                            cond_func func, void *user_data)
 { 
   if (!the_fields)
    init_fields();

   call_step_funcs(sf_list, -1.0*HUGE_VAL); // AT_BEGINNING

   while( func(the_fields, user_data) == false )
    { the_fields->step();
      call_step_funcs(sf_list, the_fields->round_time());
    };

   call_step_funcs(sf_list, +1.0*HUGE_VAL); // AT_END

 };

void meep_session::run_until(cond_func func, void *user_data)
 { run_until(step_func_list, func, user_data); }

bool cond_until(meep::fields *f, void *user_data)
{ double T = *((double *)user_data);
  return f->round_time() > T;
}

void meep_session::run_until(step_func_data *sf_list, double T)
 { run_until(sf_list, cond_until, (void *)&T); }

void meep_session::run_until(double T)
 { run_until(step_func_list, T); }

void meep_session::run_sources_plus(step_func_data *sf_list, double T)
 { if (!the_fields) init_fields();
   run_until(sf_list, the_fields->last_source_time() + T );
 }

void meep_session::run_sources_plus(double T)
 { run_sources_plus(step_func_list, T); }

void meep_session::run_sources(step_func_data *sf_list)
 { run_sources_plus(sf_list, 0.0); }

void meep_session::run_sources()
 { run_sources(step_func_list); }


typedef struct fields_decayed_data
 { 
   double dT;
   double next_check_time;
   double cur_max;
   double max_abs;
   meep::component c;
   meep::vec pt;
   double decay_by;
 } fields_decayed_data;

bool cond_fields_decayed(meep::fields *f, void *user_data)
{ 
  fields_decayed_data *data = (fields_decayed_data *)user_data;
  double absEz = abs(f->get_field(data->c, data->pt));
  data->cur_max = fmax(data->cur_max, absEz);
  if (f->round_time() <= data->next_check_time) return false;
  data->next_check_time+= data->dT;
  data->max_abs = fmax(data->max_abs, data->cur_max);
  bool status = (data->max_abs > 0.0) && (data->cur_max < data->decay_by * data->max_abs);
  data->cur_max=0.0;
  return status;
}

void meep_session::run_until_fields_decayed(step_func_data *sf_list,
                                           double dT, 
                                           meep::component c,
                                           meep::vec pt, 
                                           double decay_by)
{
  fields_decayed_data data={ dT, dT, 0.0, 0.0, c, pt, decay_by };
  run_until(sf_list, cond_fields_decayed, (void *)&data);
}

void meep_session::run_until_fields_decayed(double dT,
                                           meep::component c,
                                           meep::vec pt,
                                           double decay_by)
{ run_until_fields_decayed(step_func_list, dT, c, pt, decay_by); }

/***************************************************************/
/***************************************************************/
/***************************************************************/
bool parse_vector(char **Tokens, int nt, int NumTokens, vector3 *vec)
{
  if ( nt+4 > NumTokens ) return false;
  if (1!=sscanf(Tokens[nt+1],"%le",&(vec->x))) return false;
  if (1!=sscanf(Tokens[nt+2],"%le",&(vec->y))) return false;
  if (1!=sscanf(Tokens[nt+3],"%le",&(vec->z))) return false;

  Tokens[nt]=Tokens[nt+1]=Tokens[nt+2]=Tokens[nt+3]=0;
  return true;
}

bool parse_scalar(char **Tokens, int nt, int NumTokens, double *scalar)
{
  if ( nt+2 > NumTokens ) return false;
  if (1!=sscanf(Tokens[nt+1],"%le",scalar)) return false;
  Tokens[nt]=Tokens[nt+1]=0;
  return true;
}

/***************************************************************/
/* simple parser that accepts a string like                    */
/*  CYLINDER --center 1.0 2.0 3.0 --size 4.0 5.0 6.0           */
/* and returns a geometric_object                              */
/***************************************************************/
#define MAXSTR 1000
#define MAXTOKENS 10
geometric_object meep_session::make_object(const char *format, ...)
{
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  va_list ap;
  char description[MAXSTR];
  va_start(ap,format);
  vsnprintfEC(description,MAXSTR,format,ap);
  va_end(ap);

  char *Tokens[MAXTOKENS];
  int NumTokens=Tokenize(description, Tokens, MAXTOKENS);

  /***************************************************************/
  /* parse vector-valued and scalar-valued arguments             */
  /***************************************************************/
  vector3 center = {0.0, 0.0, 0.0};
  vector3 size   = {0.0, 0.0, 0.0};
  vector3 axis   = {0.0, 0.0, 1.0};
  vector3 e1     = {1.0, 0.0, 0.0};
  vector3 e2     = {0.0, 1.0, 0.0};
  vector3 e3     = {0.0, 0.0, 1.0};
  double radius  = 0.0;
  double radius2 = 0.0;
  double height  = HUGE_VAL;

  const char *vec_opts[]={"--center", "--size", "--e1", "--e2", "--e3", 0};
  vector3 *vecs[] ={&center,    &size,    &e1,    &e2,    &e3      };

  const char *scalar_opts[]={"--radius", "--radius2", "--height", 0};
  double *scalars[]        ={&radius,    &radius2,    &height};
 
  for(int nt=1; nt<NumTokens; nt++)
   { 
     if (Tokens[nt]==0) continue;
    
     bool handled = false;

     // try to parse vector-valued options
     for(int n=0; !handled && vec_opts[n]; n++)
      if (!strcasecmp(Tokens[nt],vec_opts[n] ))
       if (parse_vector(Tokens, nt, NumTokens, vecs[n]))
        handled=true;
       else
        meep::abort("bad values specified for option %s",Tokens[nt]);
     if (handled) continue;

     // try to parse scalar-valued options
     for(int n=0; !handled && scalar_opts[n]; n++)
      if (!strcasecmp(Tokens[nt],scalar_opts[n]))
       if (parse_scalar(Tokens, nt, NumTokens, scalars[n]))
        handled=true;
       else
        meep::abort("bad value specified for option %s",Tokens[nt]);
     if (handled) continue;

     meep::abort("unknown option %s",Tokens[nt]);
   };
     
  /***************************************************************/
  /* switch off based on the type of object requested  ***********/
  /***************************************************************/
  if ( !strcasecmp(Tokens[0], "block") )
   return make_block(meep_geom::vacuum, center, e1, e2, e3, size);
  else if ( !strcasecmp(Tokens[0], "cone") )
   return make_cone(meep_geom::vacuum, center, radius, height, axis, radius2);
  else if ( !strcasecmp(Tokens[0], "cylinder") )
   return make_cylinder(meep_geom::vacuum, center, radius, height, axis);
  else if ( !strcasecmp(Tokens[0], "ellipsoid") )
   return make_ellipsoid(meep_geom::vacuum, center, e1, e2, e3, size);
  else if ( !strcasecmp(Tokens[0], "sphere") )
   return make_sphere(meep_geom::vacuum, center, radius);
  else
   meep::abort("unknown object %s",Tokens[0]);
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
geometric_object *meep_session::add_object(material_type material,
                                           const char *format, ...)
{
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  va_list ap;
  char description[MAXSTR];
  va_start(ap,format);
  vsnprintfEC(description,MAXSTR,format,ap);
  va_end(ap);
  geometric_object o = make_object(description);
  o.material=material;
  the_object_list.items 
   = (geometric_object *)realloc( the_object_list.items,
                                 (the_object_list.num_items+1)*sizeof(o))
;
  the_object_list.items[the_object_list.num_items]=o;
  return &(the_object_list.items[the_object_list.num_items++]);
}


} // namespace meephl 
