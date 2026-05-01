/*****************************************************************************
 *
 *      Module:         RemeshRule_2B.c
 *      Description:    This module contains functions to coarsen
 *                      or refine the mesh topology.
 *
 *      Included functions:
 *              MeshCoarsen()
 *              MeshRefine()
 *              RemeshRule_2()
 *
 *****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "Home.h"
#include "Util.h"
#include "QueueOps.h"
#include "Mobility.h"
#include "ParadisSurface.h"

/*-------------------------------------------------------------------------
 *
 *      Function:    MeshCoarsen
 *      Description:
 *
 *------------------------------------------------------------------------*/
static void MeshCoarsen(Home_t *home) {
  int i, q, thisDomain, mergeDone, mergeStatus, globalOp;
  real8 tmpMag2;
  real8 cellLength, cutoffLength;
  real8 vec1x, vec1y, vec1z;
  real8 vec2x, vec2y, vec2z;
  real8 vec3x, vec3y, vec3z;
  real8 r1, r2, r3;
  real8 s, area2, areaMin, areaMin2, delta;
  real8 dvec1xdt, dvec1ydt, dvec1zdt;
  real8 dvec2xdt, dvec2ydt, dvec2zdt;
  real8 dvec3xdt, dvec3ydt, dvec3zdt;
  real8 dr1dt, dr2dt, dr3dt, dsdt, darea2dt;
  real8 newPos[3];
  real8 gp0[3], gp1[3], tmp3[3];
  Tag_t nbr1Tag, nbr2Tag, oldTag1, oldTag2, oldTag3;
  Node_t *node, *nbr, *nbr1, *nbr2, *mergedNode;
  Param_t *param;
  unsigned int iConstraintType = 0;
  Vector oConstraintVector;

  thisDomain = home->myDomain;
  param = home->param;

  cellLength = home->param->Lx / home->param->nXcells;
  cellLength = MIN(cellLength, home->param->Ly / home->param->nYcells);
  cellLength = MIN(cellLength, home->param->Lz / home->param->nZcells);

  cutoffLength = MIN(param->maxSeg, 0.45 * cellLength);

  areaMin = param->remeshAreaMin;
  areaMin2 = areaMin * areaMin;
  delta = 1.0e-16;
  /*
   *      Loop through all the nodes native to this domain looking for
   *      nodes that should be coarsened out.
   */
  for (i = 0; i < home->newNodeKeyPtr; i++) {
    node = home->nodeKeys[i];
    if (node == NULL) {
      continue;
    }
    if ((node->constraint == PINNED_NODE) ||
        (node->constraint >=SURFACE_NODE)) {
      continue;
    }
    if (node->numNbrs != 2) {
      continue;
    }
    // see if the node is inside a precipitate or not

    /*
     *          Check for various conditions that will exempt a node from
     * removal:
     *
     *          1) does not have exactly 2 arms
     *          2) node is a 'fixed' node
     *          3) node is flagged as exempt from coarsen operations
     *          4) current domain does not 'own' at least one of the  segments
     *          5) If the node's arms are on different glide planes we might
     *             not allow the node to be removed.
     */

    nbr1 = GetNeighborNode(home, node, 0);
    nbr2 = GetNeighborNode(home, node, 1);

    if ((nbr1 == NULL) || (nbr2 == NULL)) {
      printf("WARNING: Neighbor not found in MeshCoarsen\n");
      continue;
    }

    if (node->flags & NO_MESH_COARSEN) {
      continue;
    }

    if (!DomainOwnsSeg(home, OPCLASS_REMESH, thisDomain, &nbr1->myTag) &&
        !DomainOwnsSeg(home, OPCLASS_REMESH, thisDomain, &nbr2->myTag)) {
      continue;
    }

    // if the node's dynamic constraint is not a plane contraint, don't coarsen
    // it out
    ParadisSurface::GetNodeDynamicConstraint(node, iConstraintType,
                                             oConstraintVector);
    if (iConstraintType != 1) {
      continue;
    }

    /*
     *          Calculate the lengths of the node's 2 arms plus
     *          the distance between the two neighbor nodes.
     *
     *          If periodic boundaries are enabled, the nodes may
     *          be on opposite side of the problem space, so adjust
     *          the lengths/distances accordingly.
     */
    vec1x = nbr1->x - node->x;
    vec1y = nbr1->y - node->y;
    vec1z = nbr1->z - node->z;

    vec2x = nbr2->x - node->x;
    vec2y = nbr2->y - node->y;
    vec2z = nbr2->z - node->z;

    vec3x = vec2x - vec1x;
    vec3y = vec2y - vec1y;
    vec3z = vec2z - vec1z;

    GetPrimaryImage(param, &vec1x, &vec1y, &vec1z);
    GetPrimaryImage(param, &vec2x, &vec2y, &vec2z);
    GetPrimaryImage(param, &vec3x, &vec3y, &vec3z);

    r1 = sqrt(vec1x * vec1x + vec1y * vec1y + vec1z * vec1z);
    r2 = sqrt(vec2x * vec2x + vec2y * vec2y + vec2z * vec2z);
    r3 = sqrt(vec3x * vec3x + vec3y * vec3y + vec3z * vec3z);

    /*
     *          If coarsening out a node would leave a segment longer
     *          than a defined length, the node should not be removed.
     *          This 'cutoff length' is the maximum segment length
     *          if all involved nodes are within the same domain, but
     *          if a remote node is involved, we need to set the
     *          cutoff length to (at most) 1/2 the cell length.  This
     *          is needed because the remote node may potentially be involved
     *          in a simultaneous mesh coarsening in the remote domain,
     *          and although the node would not be removed, it could
     *          be repositioned resulting in a segment that spanned
     *          more than 2 cells... this is a bad thing.
     */

    // junjie
    if (home->param->EnableTwinPlaneCrossSlip == 0) {
      if (r3 > cutoffLength) {
        continue;
      }
    } else {
      if ((nbr1->csState == AT_INTERSECTION && nbr1->numNbrs == 3) ||
          (nbr2->csState == AT_INTERSECTION &&
           nbr2->numNbrs ==
               3)) { // At the intersection, it really needs coarsening
        // if (node->csState == ON_CROSS_SLIP_PLANE || node->csState ==
        // ON_SLIP_PLANE) continue;
      } else {
        if (r3 > cutoffLength) {
          continue;
        }
      }
    }
    // junjie

    /*
     *          Check if the area of the triangle defined by node
     *          and its two neighbors, plus determine if that area
     *          is increasing or decreasing.
     */

    s = 0.5 * (r1 + r2 + r3);
    area2 = (s * (s - r1) * (s - r2) * (s - r3));

    dvec1xdt = nbr1->vX - node->vX;
    dvec1ydt = nbr1->vY - node->vY;
    dvec1zdt = nbr1->vZ - node->vZ;

    dvec2xdt = nbr2->vX - node->vX;
    dvec2ydt = nbr2->vY - node->vY;
    dvec2zdt = nbr2->vZ - node->vZ;

    dvec3xdt = dvec2xdt - dvec1xdt;
    dvec3ydt = dvec2ydt - dvec1ydt;
    dvec3zdt = dvec2zdt - dvec1zdt;

    dr1dt = ((vec1x * dvec1xdt) + (vec1y * dvec1ydt) + (vec1z * dvec1zdt)) /
            (r1 + delta);
    dr2dt = ((vec2x * dvec2xdt) + (vec2y * dvec2ydt) + (vec2z * dvec2zdt)) /
            (r2 + delta);
    dr3dt = ((vec3x * dvec3xdt) + (vec3y * dvec3ydt) + (vec3z * dvec3zdt)) /
            (r3 + delta);
    dsdt = 0.5 * (dr1dt + dr2dt + dr3dt);

    darea2dt = (dsdt * (s - r1) * (s - r2) * (s - r3));
    darea2dt += s * (dsdt - dr1dt) * (s - r2) * (s - r3);
    darea2dt += s * (s - r1) * (dsdt - dr2dt) * (s - r3);
    darea2dt += s * (s - r1) * (s - r2) * (dsdt - dr3dt);

    /*
     *          If the area is less than the specified minimum and shrinking,
     *          or one of the arms is less than the minimum segment length, the
     *          node should be removed.
     */

    if (((area2 < areaMin2) && (darea2dt < 0.0)) ||
        ((r1 < param->minSeg) || (r2 < param->minSeg))) {
      mergeDone = 0;

      nbr1Tag = nbr1->myTag;
      nbr2Tag = nbr2->myTag;
      /*
       *              If either of the neighbor nodes (or any of their
       * neighbors) is in a remote domain, the operation must be treated as
       * global. This is necessary to prevent an inconsistent linkage problem.
       *              For example, given the segments A--B--C--D where nodes A,
       * B and C are in domain 1, D is in domain 2, and segment C--D is owned by
       * domain 1:  domain 1 could coarsen node B into A, then node C into D. If
       * the first operation was not communicated to the remote domain, an
       * inconsitency would arise.
       *
       *              NOTE:  It is safe to not distribute the purely local
       * coarsen operations so long as no other topological operations are done
       * after Remesh() but before the ghost node are redistributed.
       */
      globalOp = ((nbr1->myTag.domainID != thisDomain) ||
                  (nbr2->myTag.domainID != thisDomain));
      for (q = 0; q < nbr1->numNbrs; q++) {
        globalOp |= (nbr1->nbrTag[q].domainID != thisDomain);
      }

      for (q = 0; q < nbr2->numNbrs; q++) {
        globalOp |= (nbr2->nbrTag[q].domainID != thisDomain);
      }

      oldTag1 = nbr1->myTag;
      oldTag2 = node->myTag;
      oldTag3 = nbr2->myTag;
      /*
       *              If the first neighbor is not exempt from a coarsen
       *              operation, attempt to merge the nodes.
       */
      if ((nbr1->flags & NO_MESH_COARSEN) == 0) {
        newPos[X] = nbr1->x;
        newPos[Y] = nbr1->y;
        newPos[Z] = nbr1->z;
        // junjie
        /*
        if (home->param->EnableTwinPlaneCrossSlip == 1) {
            if (node->csState == AT_INTERSECTION && nbr1->csState !=
        AT_INTERSECTION) node->csState = ON_SLIP_PLANE; if (node->csState !=
        AT_INTERSECTION && nbr1->csState == AT_INTERSECTION) nbr1->csState =
        ON_SLIP_PLANE;
        }
        */
        // junjie
        mergeStatus = MergeNode(home, OPCLASS_REMESH, node, nbr1, newPos,
                                &mergedNode, globalOp);
        mergeDone = mergeStatus & MERGE_SUCCESS;
      }

      // If the merge could not be done, try using the other neighbor.
      if (mergeDone == 0) {
        if ((nbr2->flags & NO_MESH_COARSEN) == 0) {
          newPos[X] = nbr2->x;
          newPos[Y] = nbr2->y;
          newPos[Z] = nbr2->z;
          // junjie
          /*
          if (home->param->EnableTwinPlaneCrossSlip == 1) {
              if (node->csState == AT_INTERSECTION && nbr2->csState !=
          AT_INTERSECTION) node->csState = ON_SLIP_PLANE; if (node->csState !=
          AT_INTERSECTION && nbr2->csState == AT_INTERSECTION) nbr2->csState =
          ON_SLIP_PLANE;
          }
          */
          // junjie
          mergeStatus = MergeNode(home, OPCLASS_REMESH, node, nbr2, newPos,
                                  &mergedNode, globalOp);
          mergeDone = mergeStatus & MERGE_SUCCESS;
        }
      }

      // if the merge wasn't successful, continue
      if (mergeDone == 0) {
        continue;
      }

      // in some cases, the merge succeeds but then it deletes the merge node
      // after that because it is orphaned, the following check makes sure that
      // the system won't crash in this case
      if (mergedNode != NULL) {
        mergedNode->flags |= NO_MESH_COARSEN;
      }
    }
  }
}

/*-------------------------------------------------------------------------
 *
 *      Function:    TrySegBisect
 *      Description: Attempt to bisect a discretization segment during
 *                   mesh refinement by splitting the specified node.
 *                   Some restrictions apply and if any of the criteria
 *                   are not met, the segment will not be bisected.
 *
 *      Arguments:
 *          origNode  Pointer to location containing the pointer to
 *                    the node that we are considering splitting.
 *                    On return to the caller this will contain a
 *                    pointer to the node left at the position
 *                    <origNode> was in on entry to this function.
 *          nbr       Pointer to node at the end of the segment to
 *                    be bisected.
 *          armID     Index of the segment in <origNode>'s arm list
 *          segLen    Length (in units of b) of the segment to be
 *                    bisected.
 *          vec       Vector from <origNode> to <nbr>
 *          area2     area (squared) enclosed by <origNode>, <nbr> and
 *                    the other neighbor of <origNode> which is not
 *                    provided to this function.
 *          areaMax2  Maximum area (squared) needed before the segment
 *                    may be bisected (unless segLen has exceeded the
 *                    maximum permitted segment length).
 *          darea2dt  Indicates if the area indicated by <area2> is
 *                    increasing or decreasing in size.
 *          splitIsOK Additional flag indicating if the segment is
 *                    permitted to be bisected.  (Caller may have
 *                    determined the bisection should not occur).
 *                    0 == not permitted, 1 == permitted.
 *          didBisect Flag indicating if the other segment attached
 *                    to <origNode> was bisected this cycle.  Mostly for
 *                    possible future use...  0 == not bisected,
 *                    1 == bisected.
 *
 *      Returns:  1 if the segment was bisected, 0 if not
 *
 *------------------------------------------------------------------------*/
static int TrySegBisect(Home_t *home, Node_t **origNode, Node_t *nbr, int armID,
                        real8 segLen, real8 vec[3], real8 area2, real8 areaMax2,
                        real8 darea2dt, int splitIsOK, int didBisect) {
  int splitStatus, thisDomain, globalOp;
  real8 newPos[3];
  Node_t *node, *newNode;
  Param_t *param;

  thisDomain = home->myDomain;
  param = home->param;
  didBisect = 0;
  node = *origNode;

  if (!splitIsOK) {
    return (0);
  }

  /*
   *      If the segment is over the max segment length, or the area
   *      of the triangle defined by the node and its two neighbors
   *      is above the limit AND the area is increasing in size AND
   *      the segment is not too small to bisect, we cut it.
   */
  if ((segLen > param->maxSeg) ||
      ((area2 > areaMax2) && (segLen >= param->minSeg * 2.0) &&
       (darea2dt >= 0.0))) {
    /*
     *          When bisecting a segment, we always move exactly one arm
     *          from the original node to the node being created... in this
     *          case, the arm with index <armID>
     */
    newPos[X] = node->x + (vec[X] * 0.5);
    newPos[Y] = node->y + (vec[Y] * 0.5);
    newPos[Z] = node->z + (vec[Z] * 0.5);

    GetPrimaryImage(param, &newPos[X], &newPos[Y], &newPos[Z]);

    /*
     *          This should be a global operation distributed
     *          out to remote domains only if the neighbor node
     *          is in another domain.
     */
    globalOp = (nbr->myTag.domainID != node->myTag.domainID);

    newNode = NULL;
    splitStatus = SplitNode(home, node, newPos, armID, globalOp, newNode);

    if (splitStatus == SPLIT_SUCCESS) {
      didBisect = 1;
    }
  }
  return (didBisect);
}

/*-------------------------------------------------------------------------
 *
 *      Function:    MeshRefine
 *      Description:
 *
 *------------------------------------------------------------------------*/
static void MeshRefine(Home_t *home) {
  int i, thisDomain, didBisect, seg, splitIsOK;
  int splitStatus, armIndex, globalOp;
  int armID;
  int splitOK[2], splitSegList[2];
  real8 areaMax, areaMax2, maxSeg2;
  real8 delta, r1, r2, r3, s, area2, segLen;
  real8 dvec1xdt, dvec1ydt, dvec1zdt;
  real8 dvec2xdt, dvec2ydt, dvec2zdt;
  real8 dvec3xdt, dvec3ydt, dvec3zdt;
  real8 dr1dt, dr2dt, dr3dt, dsdt, darea2dt;
  real8 *vec, vec1[3], vec2[3], vec3[3];
  real8 newPos[3];
  Node_t *node, *nbr, *nbr1, *nbr2, *newNode;
  Param_t *param;

  /*
   *      Initialize some of the constants we'll need during this call.
   */
  thisDomain = home->myDomain;
  param = home->param;

  areaMax = param->remeshAreaMax;
  areaMax2 = areaMax * areaMax;
  maxSeg2 = param->maxSeg * param->maxSeg;
  delta = 1.0e-16;
  newNode = NULL;

  /*
   *      Loop through all the  native nodes looking for segments
   *      that need to be refined.
   */

  for (i = 0; i < home->newNodeKeyPtr; i++) {
    node = home->nodeKeys[i];
    if (node == NULL) {
      continue;
    }
    /*
     *          If the node has only two arms, use both arm lengths and
     *          the area in the triangle defined by the node and its
     *          neighbors to decide whether the node should be split or not.
     */
    if (node->numNbrs == 2) {
      // Calculate the lengths of the node's 2 arms plus the distance between
      // the two neighbor nodes.
      nbr1 = GetNeighborNode(home, node, 0);
      nbr2 = GetNeighborNode(home, node, 1);
      if ((nbr1 == NULL) || (nbr2 == NULL)) {
        printf("WARNING: Neighbor not found in 2-node MeshRefine\n");
        continue;
      }

      vec1[X] = nbr1->x - node->x;
      vec1[Y] = nbr1->y - node->y;
      vec1[Z] = nbr1->z - node->z;

      vec2[X] = nbr2->x - node->x;
      vec2[Y] = nbr2->y - node->y;
      vec2[Z] = nbr2->z - node->z;

      vec3[X] = vec2[X] - vec1[X];
      vec3[Y] = vec2[Y] - vec1[Y];
      vec3[Z] = vec2[Z] - vec1[Z];

      GetPrimaryImage(param, &vec1[X], &vec1[Y], &vec1[Z]);
      GetPrimaryImage(param, &vec2[X], &vec2[Y], &vec2[Z]);
      GetPrimaryImage(param, &vec3[X], &vec3[Y], &vec3[Z]);

      r1 = sqrt(vec1[X] * vec1[X] + vec1[Y] * vec1[Y] + vec1[Z] * vec1[Z]);
      r2 = sqrt(vec2[X] * vec2[X] + vec2[Y] * vec2[Y] + vec2[Z] * vec2[Z]);
      r3 = sqrt(vec3[X] * vec3[X] + vec3[Y] * vec3[Y] + vec3[Z] * vec3[Z]);

      s = 0.5 * (r1 + r2 + r3);
      area2 = (s * (s - r1) * (s - r2) * (s - r3));

      /*
       *              Determine if the area of the triangle defined by the node
       *              and its two neighbors is increasing or decreasing.
       */
      dvec1xdt = nbr1->vX - node->vX;
      dvec1ydt = nbr1->vY - node->vY;
      dvec1zdt = nbr1->vZ - node->vZ;

      dvec2xdt = nbr2->vX - node->vX;
      dvec2ydt = nbr2->vY - node->vY;
      dvec2zdt = nbr2->vZ - node->vZ;

      dvec3xdt = dvec2xdt - dvec1xdt;
      dvec3ydt = dvec2ydt - dvec1ydt;
      dvec3zdt = dvec2zdt - dvec1zdt;

      dr1dt =
          ((vec1[X] * dvec1xdt) + (vec1[Y] * dvec1ydt) + (vec1[Z] * dvec1zdt)) /
          (r1 + delta);
      dr2dt =
          ((vec2[X] * dvec2xdt) + (vec2[Y] * dvec2ydt) + (vec2[Z] * dvec2zdt)) /
          (r2 + delta);
      dr3dt =
          ((vec3[X] * dvec3xdt) + (vec3[Y] * dvec3ydt) + (vec3[Z] * dvec3zdt)) /
          (r3 + delta);
      dsdt = 0.5 * (dr1dt + dr2dt + dr3dt);

      darea2dt = (dsdt * (s - r1) * (s - r2) * (s - r3));
      darea2dt += s * (dsdt - dr1dt) * (s - r2) * (s - r3);
      darea2dt += s * (s - r1) * (dsdt - dr2dt) * (s - r3);
      darea2dt += s * (s - r1) * (s - r2) * (dsdt - dr3dt);

      /*
       *              Check if the current domain owns the segments.
       *              It may only split a segment it owns...
       */
      splitOK[0] =
          DomainOwnsSeg(home, OPCLASS_REMESH, thisDomain, &nbr1->myTag);
      splitOK[1] =
          DomainOwnsSeg(home, OPCLASS_REMESH, thisDomain, &nbr2->myTag);

      /*
       *              It would be preferable to bisect the longest segment
       *              first, so check the segment lengths and set the
       *              preferred order.
       */
      if (r1 > r2) {
        splitSegList[0] = 1;
        splitSegList[1] = 2;
      } else {
        splitSegList[0] = 2;
        splitSegList[1] = 1;
      }
      /*
       *              Loop over both segments in the preferred order and try to
       *              bisect them.
       *
       *              Note: If we bisect a segment on the first loop iteration,
       *              the following segment will only be bisected if it exceeds
       *              the maximum seg length.
       */
      didBisect = 0;
      for (seg = 0; seg < 2; seg++) {
        if (splitSegList[seg] == 1) {
          nbr = nbr1;
          vec = vec1;
          segLen = r1;
          splitIsOK = splitOK[0];
        } else {
          nbr = nbr2;
          vec = vec2;
          segLen = r2;
          splitIsOK = splitOK[1];
        }
        armID = GetArmID(node, nbr);
        didBisect = TrySegBisect(home, &node, nbr, armID, segLen, vec, area2,
                                 areaMax2, darea2dt, splitIsOK, didBisect);
      }
    } else {
      // For nodes with other than exactly two arms, bisect any arm exceeding
      // the max allowable segment length
      for (armIndex = 0; armIndex < node->numNbrs;) {
        nbr1 = GetNeighborNode(home, node, armIndex);
        if (nbr1 == NULL) {
          printf("WARNING: Neighbor not found in multinode MeshRefine\n");
          armIndex++;
          continue;
        }
        if (!DomainOwnsSeg(home, OPCLASS_REMESH, thisDomain, &nbr1->myTag)) {
          armIndex++;
          continue;
        }

        vec1[X] = nbr1->x - node->x;
        vec1[Y] = nbr1->y - node->y;
        vec1[Z] = nbr1->z - node->z;

        GetPrimaryImage(param, &vec1[X], &vec1[Y], &vec1[Z]);

        r1 = vec1[X] * vec1[X] + vec1[Y] * vec1[Y] + vec1[Z] * vec1[Z];
        if (r1 > maxSeg2) {
          newPos[X] = node->x + (vec1[X] * 0.5);
          newPos[Y] = node->y + (vec1[Y] * 0.5);
          newPos[Z] = node->z + (vec1[Z] * 0.5);

          GetPrimaryImage(param, &newPos[X], &newPos[Y], &newPos[Z]);

          /*
           *                      This should be a global operation
           *                      distributed out to remote domains only
           *                      if the segment spans domains.
           */
          globalOp = (nbr1->myTag.domainID != node->myTag.domainID);

          splitStatus =
              SplitNode(home, node, newPos, armIndex, globalOp, newNode);
        } else {
          armIndex++;
        }
      }
    }
  }
}

/*-------------------------------------------------------------------------
 *
 *      Function:    RemeshRule_2
 *      Description: Base function which invokes all subroutines
 *                   needed for handling mesh operations specific
 *                   to the second remesh rule
 *
 *------------------------------------------------------------------------*/
void RemeshRule_2(Home_t *home) {
  MeshCoarsen(home);
  MeshRefine(home);
}
