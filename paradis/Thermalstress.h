/*--------------------------------------------------------------------------
 *
 *	Thermalstress.h     Define the struct of stress
 *
 *------------------------------------------------------------------------*/

#ifndef _ThermalStress_h
#define _ThermalStress_h

struct _thermalstress {

  int n_x, n_y, n_z; /* global index in x,y,z directions */
  real8 x, y, z;     /* grid position */
  real8 s_xx, s_yy, s_zz, s_yz, s_xz, s_xy; /* stress value */
  real8 t_start, t_end;
  real8 Temp;
};

#endif
