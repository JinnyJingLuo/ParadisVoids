/*****************************************************************************
 *
 *   InData.h  Define the data structures used to hold the input data
 *
 ****************************************************************************/
#ifndef _InData_h
#define _InData_h

#include "Typedefs.h"
#include "Home.h"
#include "Tag.h"

struct _indata {
  Param_t *param;
  Thermalstress_t ****stress;
  Node_t *node; /* array of node structs */
  void *decomp; /* pointer to memory allocated to hold   */
                /* the domain decomposition read in from */
                /* the restart file(s)                   */
};

#endif
