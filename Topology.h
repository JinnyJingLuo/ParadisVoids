/***************************************************************************
 *
 *	Topology.h	Define the struct that holds all relevant data for a
 *			topological event and an event list, plus prototypes
 *			for various functions used in identifying and
 *			topological events and effecting the appropriate
 *			topology changes.
 *
 **************************************************************************/

#ifndef _TOPOLOGY_H
#define _TOPOLOGY_H

#include "Typedefs.h"
#include "Node.h"
#include "Tag.h"
#include "Vector.h"
#include "Matrix.h"
#include <stdio.h>
#include <cmath>
#include <stdlib.h>

using namespace EZ;

/*
 *      Define the set of status codes that may be returned from
 *      the MergeNode() function.
 */
#define MERGE_SUCCESS 0x01       /* Merge succeeded               */
#define MERGE_NODE_ORPHANED 0x02 /* Merge succeeded, but left an  */
                                 /* orphaned node in a remote     */
                                 /* domain                        */
#define MERGE_NO_REPOSITION 0x04 /* Merge succeeded, but the      */
                                 /* resulting node was not able   */
                                 /* to be repositioned.           */
#define MERGE_NOT_PERMITTED 0x08 /* Merge failed because the      */
                                 /* requested merge would violate */
                                 /* segment ownership rules       */
#define MERGE_DOUBLE_LINK 0x10   /* Merge failed because it would */
                                 /* have resulted in doubly       */
                                 /* connected nodes in another    */
                                 /* domain                        */
                                 /*
                                  *      Define the set of status codes that may be returned from
                                  *      the SplitNode() function.
                                  */
#define SPLIT_FAILED 0
#define SPLIT_SUCCESS 1

/*
 *      Define any processing flags that can be provided to SplitNode() to
 *      affect its behavior.
 */
#define SPLIT_DUP_SURFACE_PROP 0x01

#define OPCLASS_SEPARATION 1
#define OPCLASS_COLLISION 2
#define OPCLASS_REMESH 3
/*
 *	Prototypes for functions involved in altering topology
 */
void InitTopologyExemptions(Home_t *home);
int MergeNode(Home_t *home, int opClass, Node_t *node1, Node_t *node2,
              real8 *position, Node_t **mergedNode, int globalOp);
int NodeTopologyExemptions(Home_t *home, Node_t *node);
int RemoveDoubleLinks(Home_t *home, Node_t *node, int globalOp);
void RemoveOrphanedNodes(Home_t *home);
int SplitNode(Home_t *home, Node_t *node, real8 *position, int armIndex,
              int globalOp, Node_t *&newNode);
void DissociateBurgersVector(Home_t *home, const real8 bx, const real8 by,
                             const real8 bz, real8 *bx1, real8 *by1, real8 *bz1,
                             real8 *bx2, real8 *by2, real8 *bz2, const real8 nx,
                             const real8 ny, const real8 nz, const real8 dx,
                             const real8 dy, const real8 dz,
                             const Vector oPointToLeading);
int Dissociate(Home_t *home, Node_t *node, int nodecount, Node_t *poend1,
               Node_t *poend2, int globalOp, Node_t *&newNode);
/*int DissociateEnd(Home_t *home, Node_t *poend, const real8 bxEnd, const real8
byEnd, const real8 bzEnd, const real8 nxEnd, const real8 nyEnd, const real8
nzEnd,const real8 dx, const real8 dy, const real8 dz, const real8 vx, const
real8 vy, const real8 vz, int globalOp, Node_t*& newNode);*/
#endif /* _TOPOLOGY_H */
