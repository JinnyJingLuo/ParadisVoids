/*--------------------------------------------------------------------------
 *
 *	Node.h	Define the struct that holds all relevant data for a single
 *		node, either native or ghost
 *
 *		Notice: make depend if new files are added, otherwise always
 *			make clean after changing this file  Wei Cai 04/09/2002
 *
 *------------------------------------------------------------------------*/

#ifndef _Node_h
#define _Node_h

#include "Typedefs.h"
#include "Tag.h"

/*
 *      Define the various node 'constraints' available.  Note: these
 *      constraints are mutually exclusive.
 */
#define UNCONSTRAINED 0
#define PINNED_NODE 7
#define SURFACE_NODE 8
#define VOID_NODE 9

/*
 *      Define the bit flags that can be set for the node.  Note: these
 *      flags may be OR'ed together.
 */
#define NO_COLLISIONS 0x04
#define NO_MESH_COARSEN 0x08
#define NODE_CHK_DBL_LINK 0x10
#define NODE_NO_CROSS_SLIP 0x16

/*
 *      Used as bit flags to indicate the type of nodal data items
 *      preserved by calls to PreserveNodalData().
 */
#define NODE_POSITION 0x01
#define NODE_CURR_VEL 0x02
#define NODE_OLD_VEL 0x04

struct _node {
  int flags;

  real8 x, y, z;    /* nodal position */
  real8 fX, fY, fZ; /* nodal force: units=Pa*b^2) */
  real8 vX, vY, vZ; /* nodal velocity: units=burgMag/sec */

  real8 oldx, oldy, oldz; /* for strain increment, Moono.Rhee */
  double dNx;
  double dNy;
  double dNz;

  Tag_t myTag;

  /*
   *	nbrTag and burgID are dynamically allocated to size numNbrs, when the
   *	node is allocated.
   */
  int numNbrs;
  Tag_t *nbrTag;
  CrossSlipState csState; // qjiao: cross-slip tag (i.e., 0 ,1, 2)

  /*
   *	Arm information
   */
  real8 *armfx, *armfy, *armfz; /* arm specific force contribution */
  real8 *burgX, *burgY, *burgZ; /* burgers vector */
  real8 *nx, *ny, *nz;          /* glide plane */

  real8 *sigbLoc; /* sig.b on arms (numNbr*3) */
  real8 *sigbRem; /* sig.b on arms (numNbr*3) */
  int *piChainID;

  int constraint;

  int cellIdx;      /* cell node is currently sorted into */
  int cell2Idx;     /* cell2 node is currently sorted into */
  int cell2QentIdx; /* Index of this node's entry in the */
                    /* home->cell2QentArray.             */

  int VelocityDapmingSteps; /* 1 = native node, 0 = ghost node */

  Node_t *next; /* pointer to the next node in the queue */
                /* (ghost or free)			 */

  Node_t *nextInCell; /* used to queue node onto the current */
                      /* containing cell		       */
};

struct _nodeblock {
  NodeBlock_t *next;
  Node_t *nodes;
};

#endif
