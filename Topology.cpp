/*****************************************************************************
 *
 *	Module:		Topology.c
 *	Description:	This module contains various functions used
 *			for altering the topology of the dislocations
 *			in the problem space.
 *
 *                      In order to handle topology changes in parallel,
 *                      it is necessary to define ownership of both
 *                      nodes and segments, and adhere to some strict
 *                      rules as to what a domain may do to nodes/segments
 *                      it does own and may not do to nodes/segments it
 *                      does not own.  These rules are as follows:
 *
 *                      Ownership of a segment permits:
 *                          1) Change the node ID to which the
 *                             segment is connected
 *                          2) Change the burger's vector, force,
 *                             glide plane and all other segment
 *                             specific qualities.
 *                          3) Deletion of the segment
 *
 *                      Ownership of a node permits:
 *                          1) Change nodal position, force, velocity
 *                             and all other node specific qualities.
 *                          2) Delation of the node if (and only if)
 *                             it has no attached segments.
 *
 *	Included functions:
 *		CopyNode()
 *              DomainOwnsSeg()
 *              InitTopologyExemptions()
 *		MergeNode()
 *              RemoveDoubleLinks()
 *              RemoveOrphanedNode()
 *		SplitNode()
 *
 *****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include "Home.h"
#include "Util.h"
#include "Comm.h"
#include "Mobility.h"
#include "Matrix.h"
#include "Vector.h"

typedef struct {
  real8 f1x, f2x;
  real8 f1y, f2y;
  real8 f1z, f2z;
} SegData_t;

/*---------------------------------------------------------------------------
 *
 *	Function:	RemoveOrphanedNodes
 *	Description:	This function scans the array of nodes local to
 *                      this domain looking for any nodes that no longer
 *                      have any arms.   And such nodes are removed
 *                      locally.
 *
 *                      Note: The remote domains are NOT notified of this
 *                      node removal at this time.  This should not be
 *                      a problem since the ghost node communication
 *                      will essentially update the remote domains with
 *                      the proper node list.
 *
 *-------------------------------------------------------------------------*/
void RemoveOrphanedNodes(Home_t *home) {
  int i;
  Node_t *node;
  for (i = 0; i < home->newNodeKeyPtr; i++) {
    node = home->nodeKeys[i];
    if (node == (Node_t *)NULL)
      continue;
    if (node->numNbrs == 0) {
      RemoveNode(home, node, 0);
    }
  }
}

/*---------------------------------------------------------------------------
 *
 *      Function:       InitTopologyExemptions
 *      Description:    For all nodes (including ghost nodes) clear all
 *                      flags that would exempt a node and its arms from
 *                      any topological changes.
 *
 *-------------------------------------------------------------------------*/
void InitTopologyExemptions(Home_t *home) {
  int i;
  Node_t *node;
  for (i = 0; i < home->newNodeKeyPtr; i++) {
    node = home->nodeKeys[i];
    if (node == NULL) {
      continue;
    }
    node->flags &= ~(NO_COLLISIONS | NO_MESH_COARSEN);
  }
  node = home->ghostNodeQ;
  while (node != NULL) {
    node->flags &= ~(NO_COLLISIONS | NO_MESH_COARSEN);
    node = node->next;
  }
}

/*---------------------------------------------------------------------------
 *
 *      Function:     RemoveDoubleLinks
 *      Description:  This function looks for double links between
 *                    the specified node and any of its neighbors.
 *                    If a double link is found, the links are either
 *                    removed (if the burgers vectors cancel) or merged
 *                    into a single link.  If any actions performed
 *                    by this function result in an orphaned neighbor
 *                    node, the neighbor will be removed iff the neighbor
 *                    is local to the current domain.
 *
 *      Arguments:
 *          node      pointer to local node to be examined.
 *          globalOp  Flag indicating if this is a global operation
 *                    that should be added to the list of ops
 *                    distributed to neighboring domains.
 *
 *      Returns:    1 if a node was orphaned by an operation
 *                  0 in all other cases
 *
 *-------------------------------------------------------------------------*/
int RemoveDoubleLinks(Home_t *home, Node_t *node, int globalOp) {
  int i, j, nbrIndex, nbrDomain, thisDomain;
  int arm, needGlide;
  int nodeOrphaned = 0;
  real8 bx, by, bz;
  real8 nx, ny, nz;
  double dX = 0.0;
  double dY = 0.0;
  double dZ = 0.0;
  real8 bsum[3];
  Node_t *nbrNode;
  real8 burg[3], glidePlane[3];

  thisDomain = home->myDomain;

  for (i = 0; i < (node->numNbrs - 1); i++) {
    nbrDomain = node->nbrTag[i].domainID;
    nbrIndex = node->nbrTag[i].index;
    for (j = i + 1; j < node->numNbrs; j++) {
      if ((node->nbrTag[j].domainID == nbrDomain) &&
          (node->nbrTag[j].index == nbrIndex)) {
        /*
         *                  Found redundant arms linking 2 nodes.  If
         *                  the burgers vectors cancel, remove both links
         *                  otherwise, reset the burgers vector for the
         *                  first arm and remove the second.
         */
        bx = node->burgX[i] + node->burgX[j];
        by = node->burgY[i] + node->burgY[j];
        bz = node->burgZ[i] + node->burgZ[j];

        /*
         *                  Just a safety check to prevent tiny non-zero
         *                  components of burgers vector due to machine
         *                  precision or round-off issues.
         */
        if (fabs(bx) < 1.0e-06)
          bx = 0.0;
        if (fabs(by) < 1.0e-06)
          by = 0.0;
        if (fabs(bz) < 1.0e-06)
          bz = 0.0;

        if (fabs((bx * bx) + (by * by) + (bz * bz)) < 1.0e-6) {
          nbrNode = GetNeighborNode(home, node, j);

          if (nbrNode == NULL) {
            Fatal("Neighbor not found at %s line %d\n", __FILE__, __LINE__);
          }

          ChangeArmBurg(home, node, &nbrNode->myTag, 0.0, 0.0, 0.0, 0.0, 0.0,
                        0.0, globalOp, DEL_SEG_HALF, true);
          ChangeArmBurg(home, nbrNode, &node->myTag, 0.0, 0.0, 0.0, 0.0, 0.0,
                        0.0, globalOp, DEL_SEG_HALF, true);
          /*
           *                      It's possible the neighbor node no longer has
           * any arms.  If the neighbor is local, go ahead and remove it,
           * otherwise, just set a status flag letting the caller know a remote
           * node has been orphaned.
           */
          if (nbrNode->numNbrs == 0) {
            if (nbrNode->myTag.domainID == thisDomain) {
              RemoveNode(home, nbrNode, globalOp);
            } else {
              nodeOrphaned = 1;
            }
          }

          /*
           *                      We removed multiple links, so set i back by
           * one so we don't skip over an arm in the outer loop.
           */
          i--;
        } else {
          // the slip plane normal for the new segment should be normal to the
          // burgers vector and the line direction between the two nodes
          nbrNode = GetNeighborNode(home, node, j);
          if (nbrNode == NULL) {
            Fatal("Neighbor not found at %s line %d\n", __FILE__, __LINE__);
          }
          dX = nbrNode->x - node->x;
          dY = nbrNode->y - node->y;
          dZ = nbrNode->z - node->z;
          GetPrimaryImage(home->param, &dX, &dY, &dZ);
          Vector oDirection(dX, dY, dZ);
          Vector oBurgers(bx, by, bz);
          Vector oNormal = oDirection ^ oBurgers;
          // junjie
          real8 segPlaneX, segPlaneY, segPlaneZ;
          real8 A, B, C, len;
          A = home->param->A;
          B = home->param->B;
          C = home->param->C;
          len = sqrt(pow(A, 2.0) + pow(B, 2.0) + pow(C, 2.0));
          segPlaneX = A / len;
          segPlaneY = B / len;
          segPlaneZ = C / len;
          Vector InPutNormal(segPlaneX, segPlaneY, segPlaneZ);
          double dTolerance = 1.0E-6;
          // junjie

          if (oNormal.Length() < dTolerance) {
            // the cross product is zero, the segment is a screw dislocation
            burg[0] = bx;
            burg[1] = by;
            burg[2] = bz;
            PickScrewGlidePlane(home, burg, glidePlane);
            nx = glidePlane[0];
            ny = glidePlane[1];
            nz = glidePlane[2];
          } else {
            // junjie
            /*
             *                       Just a safety check to prevent tiny
             * non-zero components of burgers vector due to machine precision or
             * round-off issues.
             */
            if (home->param->EnableTwinPlaneCrossSlip == 1 &&
                (oNormal ^ InPutNormal).Length() < dTolerance) {
              oNormal = InPutNormal;
            }

            if (fabs(oNormal.GetX()) < 1.0e-06)
              oNormal.SetX(0.0);
            if (fabs(oNormal.GetY()) < 1.0e-06)
              oNormal.SetY(0.0);
            if (fabs(oNormal.GetZ()) < 1.0e-06)
              oNormal.SetZ(0.0);
            // junjie
            oNormal.Normalize();
            nx = oNormal.GetX();
            ny = oNormal.GetY();
            nz = oNormal.GetZ();
          }
          // since we are merging two segments into one nonzero segments, we
          // need to properly handle their chain IDs. If the two segments belong
          // to the same chain, the chain ID of the merged segment will be the
          // common segments ID, otherwise, it will be zero for now
          int iChainID = 0;
          if (node->piChainID[i] == node->piChainID[j]) {
            iChainID = node->piChainID[i];
          }
          ChangeArmBurg(home, node, &nbrNode->myTag, bx, by, bz, nx, ny, nz,
                        globalOp, DEL_SEG_HALF, false);
          ChangeArmBurg(home, nbrNode, &node->myTag, -bx, -by, -bz, nx, ny, nz,
                        globalOp, DEL_SEG_HALF, false);
          SetChainID(node, nbrNode, iChainID);
        }
        /*
         *                  Redundant links removed (or merged), so break
         *                  out of the inner loop and continue the outer loop.
         */
        break;

      } /* if redundant links found */
    }
  }

  return (nodeOrphaned);
}
/*---------------------------------------------------------------------------
 *
 *	Function:	SplitNode
 *	Description:	Create a new node and transfer the specified set of
 *			connections from an existing to to the new node.  If
 *			necessary, a new link will also be created between the
 *			existing node and the new node.
 *
 *	Arguments:
 *		node            Pointer to the node to be split
 *              pos1            coordinates at which splitNode1 will be
 *                              left after the split
 *              pos2            coordinates at which splitNode2 will be
 *                              left after the split
 *              vel1            velocity assigned to splitNode1 after the split
 *              vel2            velocity assigned to splitNode2 after the split
 *		armCount	number of arms of the original node selected
 *                              to be split off.
 *		armList		pointer to array of integers indicating the
 *				arms of existing node that are to be split off
 *		globalOp	Flag indicating if this is a global operation
 *				that should be added to the list of ops
 *				distributed to neighboring domains.
 *              splitNode1      ptr to ptr to node to which all unselected arms
 *                              of the original node will be attached after the
 *                              split.  Returned to caller.
 *              splitNode2      ptr to ptr to node to which all selected arms
 *                              of the original node will be attached after the
 *                              after the split.  Returned to caller.
 *
 *	Returns:		1 if the split was successful
 *                              0 in all other cases
 *
 *-------------------------------------------------------------------------*/
int SplitNode(Home_t *home, Node_t *node, real8 *position, int armIndex,
              int globalOp, Node_t *&newNode) {
  int i;
  real8 bx, by, bz;
  real8 nx, ny, nz;
  real8 segPlaneX, segPlaneY, segPlaneZ;
  real8 eps = 1.0e-12;
  real8 ftmp[3];
  real8 burg[3];
  Node_t *nbrNode;
  Param_t *param;
  int thisDomain;

  thisDomain = home->myDomain;
  param = home->param;

  // basic checks
  if (node->myTag.domainID != thisDomain) {
    return (SPLIT_FAILED);
  }

  /*
   *      If we are only moving a single segment from the original node to the
   *      new node, it (likely) means we're simply bisecting a segment during
   *      MeshRefine().  In that case, preserve the segment's glide plane
   *      so when we add a new segment between the original node and new node
   *      we can have the new segment properly inherit the glide plane.
   */
  segPlaneX = node->nx[armIndex];
  segPlaneY = node->ny[armIndex];
  segPlaneZ = node->nz[armIndex];

  /*
   *	Add a new node.  Warning:  GetNewNativeNode() pops a node off
   *      the queue, but does not free the arm related arrays.  It does
   *      this on the assumption that the number of arms previously
   *      assigned to the node is likely to be the same that will
   *      be reattached to the node and that the code will overwrite
   *      the existing arm data.  The following code, however, assumes
   *	the new node has zero arms initially, so we MUST free those arrays
   *	and reset the arm count!
   */
  newNode = GetNewNativeNode(home);
  FreeNodeArms(newNode);

  newNode->constraint = UNCONSTRAINED;

  newNode->dNx = 0.0;
  newNode->dNy = 0.0;
  newNode->dNz = 0.0;

  // set the new node position, old position and velocity
  newNode->x = position[X];
  newNode->y = position[Y];
  newNode->z = position[Z];

  newNode->oldx = newNode->x;
  newNode->oldy = newNode->y;
  newNode->oldz = newNode->z;

  /*
   *      If this is a global operation that remote domains will have to
   *      know about, add the operation to the list that will be sent
   *      to the neighboring domains.  NOTE: Certain split operations
   *      may have altered the coordinates of <origNode>, so also send
   *      out updated coordinates to the remote domains.
   */
  if (globalOp) {
    AddOp(home, SPLIT_NODE, node->myTag.domainID, node->myTag.index,
          newNode->myTag.domainID, newNode->myTag.index, -1,
          -1,            /* 3rd node tag unneeded */
          0.0, 0.0, 0.0, /* burgers vector unneeded*/
          newNode->x, newNode->y, newNode->z, 0.0, 0.0,
          0.0); // normal not needed
    /* arguments for velocity data*/
    AddOp(home, RESET_COORD, node->myTag.domainID, node->myTag.index, -1,
          -1,            /* 2nd tag unneeded */
          -1, -1,        /* 3rd node tag unneeded */
          0.0, 0.0, 0.0, /* burgers vector unneeded */
          node->x, node->y, node->z, 0.0, 0.0,
          0.0); /* plane normal not needed */
  }

  /*
   *	For each arm of the original node in armList, get the neighbor node,
   *	change the nbrNode->origNode linkage to nbrNode->newNode, and
   *	remove the origNode->nbrNode linkage.  Note: the DEL_SEG_NONE flag
   *	is passed to ChangeArmBurg() to prevent the function from treating
   *	the segments as annihilated since they are simply having an endpoint
   *	shifted.
   */
  nbrNode = GetNeighborNode(home, node, armIndex);
  if (nbrNode == (Node_t *)NULL) {
    Fatal("Neighbor not found at %s line %d\n", __FILE__, __LINE__);
  }

  ftmp[X] = node->armfx[armIndex];
  ftmp[Y] = node->armfy[armIndex];
  ftmp[Z] = node->armfz[armIndex];

  // save the chain ID before removing the arm
  int iChainID = node->piChainID[armIndex];
  SubtractSegForce(home, node, nbrNode);
  GetBurgersVectorNormal(home, node, nbrNode, &bx, &by, &bz, &nx, &ny, &nz);
  ChangeConnection(home, nbrNode, &node->myTag, &newNode->myTag, globalOp);
  InsertArm(home, newNode, &node->nbrTag[armIndex], bx, by, bz, nx, ny, nz,
            iChainID, globalOp);
  ChangeArmBurg(home, node, &nbrNode->myTag, 0, 0, 0, 0, 0, 0, globalOp,
                DEL_SEG_NONE, true);
  ResetSegForces(home, newNode, &nbrNode->myTag, ftmp[X], ftmp[Y], ftmp[Z],
                 globalOp);

  /*
   *	We will probably be re-evaluating forces on this new node.  In
   *      order to do that, the new node must be assigned to the proper
   *      cell.
   */
  AssignNodeToCell(home, newNode);

  /*
   *	If necessary, create a link between the original node and new node
   *	in order to conserve burgers vector.  First sum up the burgers vectors
   *	for the arms of the new node...
   *
   */
  bx = 0.0;
  by = 0.0;
  bz = 0.0;
  for (i = 0; i < newNode->numNbrs; i++) {
    bx += newNode->burgX[i];
    by += newNode->burgY[i];
    bz += newNode->burgZ[i];
  }

  /*
   *      Under two circumstances we need to add an additional link between
   *      the two nodes:
   *      1) if the sum of burgers vectors on the new node is non-zero
   *      2) if the sum of burgers vectors on the new node is zero BUT
   *         the node is 'pinned';  Scenario is this:
   *         -  the single arm of a pinned node is being split
   *         -  the 'pinned' node is repositioned and unpinned
   *         -  the new node is positioned at the location of
   *            the original node and pinned.
   *         -  the original node maintains its single connection
   *
   *         At this point, the new node is pinned, but has no arms
   *         hence the sum of its burgers vectors is zero, so no
   *         link would be created between the original node and the
   *         new node.  This results in the pinned node being orphaned
   *         and deleted, and the other node is no longer attached to
   *         a pinned node and this is bad.
   */
  if (fabs(bx * bx + by * by + bz * bz) > 0.001) {
    // use the original glide plane normal when constructing the new segment
    InsertArm(home, node, &newNode->myTag, bx, by, bz, segPlaneX, segPlaneY,
              segPlaneZ, iChainID, globalOp);
    InsertArm(home, newNode, &node->myTag, -bx, -by, -bz, segPlaneX, segPlaneY,
              segPlaneZ, iChainID, globalOp);
  }
  return (SPLIT_SUCCESS);
}

/*---------------------------------------------------------------------------
 *
 *	Function:	MergeNode
 *	Description:	This function 'merges' two nodes by moving all
 *			arms of <deadNode> to <targetNode>, and then
 *			completely removing <deadNode>.  If the merge
 *			resulted in self-links from <targetNode> back
 *			to itself, or multiple links from <targetNode>
 *			to any other single node, these redundant links
 *			will be dealt with before returning to the caller.
 *
 *	Arguments:
 *              opClass         Flag indicating the class of topological
 *                              operation invoking this function.
 *                              initiated.  Valid types are:
 *
 *                                  OPCLASS_SEPARATION
 *                                  OPCLASS_COLLISION
 *                                  OPCLASS_REMESH
 *
 *              node1           pointer to first node to be merged
 *              node2           pointer to second node to be merged
 *              position        coordinates (x,y,z) at which final merged
 *                              node is to be placed.
 *              mergedNode      pointer to location in which to return
 *                              pointer to the node resulting from the merge.
 *                              A NULL pointer will be returned if the
 *                              the merge fails.
 *              status          pointer to location in which to return
 *                              a completion status to the caller.  Valid
 *                              statuses are the following, where all statuses
 *                              indicating success may be logically OR'ed
 *
 *                                  MERGE_SUCCESS
 *                                  MERGE_NO_REPOSITION
 *                                  MERGE_NODE_ORPHANED
 *                                  MERGE_NOT_PERMITTED
 *                                  MERGE_DOUBLE_LINK
 *
 *		globalOp	Flag indicating if this is a global operation
 *				that should be added to the list of ops
 *				distributed to neighboring domains.
 *
 *-------------------------------------------------------------------------*/
int MergeNode(Home_t *home, int opClass, Node_t *node1, Node_t *node2,
              real8 *position, Node_t **mergedNode, int globalOp) {
  int i;
  int armID, armID1, armID2, thisDomain;
  int node1Deletable, node2Deletable;
  int targetIsLocal, nbrIsLocal, targetOwnsSeg;
  real8 bx, by, bz;
  real8 nx, ny, nz;
  real8 ftmp[3];
  Node_t *nbrNode, *targetNode, *deadNode;
  Tag_t tmpTag, *nbr1Tag, *nbr2Tag;

  thisDomain = home->myDomain;

  targetNode = node1;
  deadNode = node1;
  *mergedNode = (Node_t *)NULL;
  int iStatus = MERGE_NOT_PERMITTED;

  /*
   *      Check if we are permitted to delete node1.  In order
   *      for a node to be deletable, it must satisfy the following:
   *
   *      1) Node may not be a 'fixed' node
   *      2) Node must be owned by the current domain
   *      3) All arms of the node must be owned by the current domain
   *      4) Node may not be exempt from (additional) collisions this cycle
   *      5) Node may not be exempt from (additional) remeshes this cycle
   *      6) If node is a surface node, *both* nodes must be on the surface
   *         for it to be deletable.  (May need to add restrictions
   *         that surface nodes be on the same surface)
   */
  node1Deletable = (node1->myTag.domainID == thisDomain);
  node1Deletable &= (node1->constraint != PINNED_NODE);
  node1Deletable &= (node1->constraint != SURFACE_NODE);
  node1Deletable &= ((node1->flags & NO_COLLISIONS) == 0);
  node1Deletable &= ((node1->flags & NO_MESH_COARSEN) == 0);
  // junjie
  // if (home->param->EnableTwinPlaneCrossSlip == 1) {
  //    if (node1->csState == AT_INTERSECTION && node2->csState ==
  //    AT_INTERSECTION && node1->numNbrs == 3 && node2->numNbrs == 3) ; else
  //    node1Deletable &= (node1->csState != AT_INTERSECTION);
  //}
  // junjie

  for (armID = 0; armID < node1->numNbrs; armID++) {
    node1Deletable &=
        DomainOwnsSeg(home, opClass, thisDomain, &node1->nbrTag[armID]);
  }

  /*
   *      If we can delete node1, use node2 as the target, but if node1
   *      cannot be removed, use node2 if we can. If node2 cannot
   *      be deleted either, return a failure to the caller.
   */
  if (node1Deletable) {
    targetNode = node2;
    deadNode = node1;
  } else {
    node2Deletable = (node2->myTag.domainID == thisDomain);
    node2Deletable &= (node2->constraint != PINNED_NODE);
    node2Deletable &= (node2->constraint != SURFACE_NODE);
    node2Deletable &= ((node2->flags & NO_COLLISIONS) == 0);
    node2Deletable &= ((node2->flags & NO_MESH_COARSEN) == 0);
    // junjie
    // if (home->param->EnableTwinPlaneCrossSlip == 1) {
    //    if (node1->csState == AT_INTERSECTION && node2->csState ==
    //    AT_INTERSECTION && node1->numNbrs == 3 && node2->numNbrs == 3) ; else
    //    node2Deletable &= (node2->csState != AT_INTERSECTION);
    //}
    // junjie

    for (armID = 0; armID < node2->numNbrs; armID++) {
      node2Deletable &=
          DomainOwnsSeg(home, opClass, thisDomain, &node2->nbrTag[armID]);
    }

    if (!node2Deletable) {
      *mergedNode = NULL;
      return MERGE_NOT_PERMITTED;
    }
    targetNode = node1;
    deadNode = node2;
  }
  // junjie
  /*
  if (home->param->EnableTwinPlaneCrossSlip == 1) {
      if (node1->csState == AT_INTERSECTION || node2->csState ==
  AT_INTERSECTION) { node1->csState = ON_SLIP_PLANE; node2->csState =
  ON_SLIP_PLANE;
      }
  }
  */

  /*
   *      Coarsening a node out *may* leave a double connection.  That
   *      is only permitted if the current domain would own both
   *      links in the resulting double link and hence be able to
   *      collapse/annihilate the links.
   */
  targetIsLocal = (targetNode->myTag.domainID == thisDomain);
  for (armID1 = 0; armID1 < targetNode->numNbrs; armID1++) {
    nbr1Tag = &targetNode->nbrTag[armID1];
    nbrIsLocal = (nbr1Tag->domainID == thisDomain);
    for (armID2 = 0; armID2 < deadNode->numNbrs; armID2++) {
      nbr2Tag = &deadNode->nbrTag[armID2];
      /*
       *              If the target node and dead node both have links to a
       *              3rd node, check the ownership of the link between the
       *              target node and the 3rd node...
       */
      if ((nbr1Tag->domainID == nbr2Tag->domainID) &&
          (nbr1Tag->index == nbr2Tag->index)) {
        if (targetIsLocal && nbrIsLocal) {
          continue;
        }

        targetOwnsSeg =
            DomainOwnsSeg(home, opClass, targetNode->myTag.domainID, nbr1Tag);

        if ((targetIsLocal && targetOwnsSeg) ||
            (nbrIsLocal && !targetOwnsSeg)) {
          continue;
        }

        /*
         *                  Link is remotely owned, so we cannot risk creating
         *                  a double link here...
         */
        *mergedNode = NULL;
        return MERGE_DOUBLE_LINK;
      }
    }
  }

  /*
   *      We've made it past all the checks that would have prevented the
   *      merge operation from being performed.  The target and dead nodes
   *      have been selected, so if the target node is local to this
   *      domain, go ahead and reposition it.  If we can't reposition the
   *      node, the merge fails
   */

  // at this point, store the coordinates of the merging nodes and the target
  // position in order to use them later to check for precipitate shearing we
  // need to make 2 checks, one for the triangle neighbour - dead - old target
  // and the other for the triangle neighbour - dead - new target because
  // changing the target location sweeps some area as well
  Polyhedron *poContainingPrecipitate = NULL;
  Polyhedron *poTempPrecipitate = NULL;
  Point oDeadPoint;
  Point oOldTargetPoint;
  Point oNewTargetPoint;
  Point oNeighbourPoint;
  Vector oBurgersVector;

  bool bIsInPrecipitate = home->poPrecipitateServer->IsSegmentShearCapable(
      deadNode, targetNode, poContainingPrecipitate);
  oDeadPoint.Set(deadNode->x, deadNode->y, deadNode->z);
  oOldTargetPoint.Set(targetNode->x, targetNode->y, targetNode->z);
  // the third point and the Burgers vector will be set later

  bool bIsNodeMobile = (targetNode->constraint != PINNED_NODE) &&
                       (targetNode->constraint != SURFACE_NODE);
  if ((targetNode->myTag.domainID == thisDomain) && bIsNodeMobile) {
    targetNode->x = position[X];
    targetNode->y = position[Y];
    targetNode->z = position[Z];
    oNewTargetPoint.Set(targetNode->x, targetNode->y, targetNode->z);
    // update the old coordinates of the target node
    targetNode->oldx = targetNode->x;
    targetNode->oldy = targetNode->y;
    targetNode->oldz = targetNode->z;
    /*
     *          If we're updating the location of the target node, we may
     *          have to communicate the new position to remote domains.
     */
    if (globalOp) {
      AddOp(home, RESET_COORD, targetNode->myTag.domainID,
            targetNode->myTag.index, -1, -1, /* 2nd tag unneeded */
            -1, -1,                          /* 3rd node tag unneeded */
            0.0, 0.0, 0.0,                   /* burgers vector unneeded */
            targetNode->x, targetNode->y, targetNode->z, 0.0, 0.0,
            0.0); /* plane normal not needed */
    }
    *mergedNode = targetNode;
    iStatus = MERGE_SUCCESS;
  } else {
    *mergedNode = NULL;
    return MERGE_NO_REPOSITION;
  }

  /*
   *      If there are any connections from the target/destination node
   *      to the dead node, use ChangeArmBurg() to get rid of those
   *      connections, then move all connections from the dead node to
   *      the target node, and add a new connection from the target node
   *      to each of the dead node's neighbors.
   */
  ChangeArmBurg(home, targetNode, &deadNode->myTag, 0.0, 0.0, 0.0, 0.0, 0.0,
                0.0, globalOp, DEL_SEG_NONE, true);
  ChangeArmBurg(home, deadNode, &targetNode->myTag, 0.0, 0.0, 0.0, 0.0, 0.0,
                0.0, globalOp, DEL_SEG_NONE, true);

  /*
   *      Move all connections from the dead node to the target node
   *      and add a new connection from the target node to each of the
   *      dead node's neighbors.
   */
  // also check for precipitate shearing at this step
  bool bIsInPrecipitateLocal = false;
  Vector oNormal;
  int iChainID = 0;
  for (i = (deadNode->numNbrs - 1); i >= 0; i--) {
    tmpTag.domainID = deadNode->nbrTag[i].domainID;
    tmpTag.index = deadNode->nbrTag[i].index;
    nbrNode = GetNeighborNode(home, deadNode, i);

    bIsInPrecipitateLocal = bIsInPrecipitate ||
                            home->poPrecipitateServer->IsSegmentShearCapable(
                                deadNode, nbrNode, poTempPrecipitate);
    oNeighbourPoint.Set(nbrNode->x, nbrNode->y, nbrNode->z);
    oNeighbourPoint = GetNearestImage(home->param, deadNode, oNeighbourPoint);
    oBurgersVector.Set(deadNode->burgX[i], deadNode->burgY[i],
                       deadNode->burgZ[i]);
    oNormal.Set(deadNode->nx[i], deadNode->ny[i], deadNode->nz[i]);

    if (nbrNode == NULL) {
      Fatal("Neighbor not found at %s line %d\n", __FILE__, __LINE__);
    }

    ftmp[X] = deadNode->armfx[i];
    ftmp[Y] = deadNode->armfy[i];
    ftmp[Z] = deadNode->armfz[i];

    // save the chain ID before removing the dead arm
    iChainID = deadNode->piChainID[i];
    GetBurgersVectorNormal(home, deadNode, nbrNode, &bx, &by, &bz, &nx, &ny,
                           &nz);
    ChangeConnection(home, nbrNode, &deadNode->myTag, &targetNode->myTag,
                     globalOp);
    ChangeArmBurg(home, deadNode, &nbrNode->myTag, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                  globalOp, DEL_SEG_NONE, true);
    InsertArm(home, targetNode, &tmpTag, bx, by, bz, nx, ny, nz, iChainID,
              globalOp);
    ResetSegForces(home, targetNode, &tmpTag, ftmp[X], ftmp[Y], ftmp[Z],
                   globalOp);

    StoreMergePlasticSweep(home, oNeighbourPoint, oDeadPoint, oOldTargetPoint,
                           oBurgersVector, oNormal);
    // now we have all what we need to check for precipitate shearing
    if (bIsInPrecipitateLocal) {
      home->poPrecipitateServer->UpdateTriangleSweep(
          home, oNeighbourPoint, oDeadPoint, oOldTargetPoint, oBurgersVector);
    }
  }

  // if the merge happens inside a precipitate, we need to check for the target
  // neighbour sweeps as well
  bIsInPrecipitateLocal = false;
  if (bIsInPrecipitate) {
    for (i = 0; i < targetNode->numNbrs; i++) {
      nbrNode = GetNeighborNode(home, targetNode, i);

      bIsInPrecipitateLocal = bIsInPrecipitate ||
                              home->poPrecipitateServer->IsSegmentShearCapable(
                                  targetNode, nbrNode, poTempPrecipitate);
      oNeighbourPoint.Set(nbrNode->x, nbrNode->y, nbrNode->z);
      oNeighbourPoint =
          GetNearestImage(home->param, targetNode, oNeighbourPoint);
      oBurgersVector.Set(targetNode->burgX[i], targetNode->burgY[i],
                         targetNode->burgZ[i]);
      oNormal.Set(targetNode->nx[i], targetNode->ny[i], targetNode->nz[i]);
      StoreMergePlasticSweep(home, oNeighbourPoint, oOldTargetPoint,
                             oNewTargetPoint, oBurgersVector, oNormal);
      // now we have all what we need to check for precipitate shearing
      if (bIsInPrecipitateLocal) {
        home->poPrecipitateServer->UpdateTriangleSweep(
            home, oNeighbourPoint, oOldTargetPoint, oNewTargetPoint,
            oBurgersVector);
      }
    }
  }

  /*
   *      The <deadNode> should have no more arms, so we can safely
   *      remove it now.
   */
  RemoveNode(home, deadNode, globalOp);
  /*
   *      Check for multiple links between the target node and any other
   *      single node.  If found, the burgers vectors for the multiple links
   *      will either cancel each other out in which case both links are
   *      removed, or the links get combined into a single link.
   */

  if (RemoveDoubleLinks(home, targetNode, globalOp)) {
    iStatus |= MERGE_NODE_ORPHANED;
  }

  /*
   *      It's possible that the target node no longer has any arms.  If
   *      it is local, go ahead and remove it, otherwise, just
   *      set a status flag letting the caller know a remote node has been
   *      orphaned.
   */
  if (targetNode->numNbrs == 0) {
    if (targetNode->myTag.domainID == thisDomain) {
      RemoveNode(home, targetNode, globalOp);
      *mergedNode = NULL;
    } else {
      iStatus |= MERGE_NODE_ORPHANED;
    }
  }

  return iStatus;
}

/*---------------------------------------------------------------------------
 *
 *	Function:	Dissociate
 *	Description:	This function is used in the cross-slip on the twin
 *plane, FCC dislocations are approximate to be perfect in bulk. However, on the
 *twin boundary, the stacking fault engs are low( Jaafar ). So partials tend to
 *push each other away. perfect screws will dissociate on twin plane. This
 *function will dissociate a screw segment that is going to cross slip onto the
 *twin plane, the input node should only have 2(surface node(todo): just
 *dissociate the node; bulk nodes: both of the arms need to dissociate) or 3
 *neighbors(1 of the 3 arms need to dissociate)(junjie: need to be careful about
 *the case when 3-arm nodes hitting the TP) Junjie: The forces on the input
 *nodes are not handled properly, if you are using node forces after this
 *founction, can potentially cause a problem.
 *
 *	Arguments:
 *		node            Pointer to the node to dissociate
 *              pos1 leading           coordinates at which dissociateNode1 will
 *be left after the dissociate pos2 trailing           coordinates at which
 *dissociateNode2 will be left after the dissociate armCount	number of arms
 *of the original node selected to be dissociated. globalOp	Flag indicating
 *if this is a global operation that should be added to the list of ops
 *				distributed to neighboring domains.
 *      newNode     ptr point to the new node from the dissociation
 *
 *	Returns:		1 if the dissociate was successful
 *                              0 in all other cases
 *
 *-------------------------------------------------------------------------*/
int Dissociate(Home_t *home, Node_t *node, int nodecount, Node_t *poend1,
               Node_t *poend2, int globalOp, Node_t *&newNode) {
  int freq = 200;
  int iChainID = 0;
  int i;
  real8 bx, by, bz;
  real8 bx1, by1, bz1;
  real8 bx2, by2, bz2;
  real8 dx, dy, dz;
  real8 nx, ny, nz;
  real8 segPlaneX, segPlaneY, segPlaneZ;
  real8 pos[3];
  real8 A, B, C, len;
  real8 lenDissociate = home->param->rmax + 1.0;
  Node_t *nbrNode1, *nbrNode2, *curNode, *curNode1, *curNode2, *lastNode;
  real8 dTol = 1E-6;
  A = home->param->A;
  B = home->param->B;
  C = home->param->C;
  len = sqrt(pow(A, 2.0) + pow(B, 2.0) + pow(C, 2.0));
  segPlaneX = A / len;
  segPlaneY = B / len;
  segPlaneZ = C / len;

  // basic checks: only when this entire segment is inside this domain, can it
  // dissociate
  if (poend1->myTag.domainID != home->myDomain ||
      poend2->myTag.domainID != home->myDomain) {
    return (SPLIT_FAILED);
  }

  nbrNode1 = GetNeighborNode(home, node, 0);
  nbrNode2 = GetNeighborNode(home, node, 1);
  if (nbrNode2 != poend2) {
    nbrNode1 = nbrNode2;
    nbrNode2 = poend2;
  }
  if (nbrNode1 == (Node_t *)NULL || nbrNode2 == (Node_t *)NULL) {
    Fatal("Neighbor not found at %s line %d\n", __FILE__, __LINE__);
  }
  GetBurgersVectorNormal(home, node, nbrNode2, &bx, &by, &bz, &nx, &ny, &nz);
  dx = nbrNode2->x - node->x;
  dy = nbrNode2->y - node->y;
  dz = nbrNode2->z - node->z;
  Vector ovelo(node->vX, node->vY, node->vZ);
  Vector oNormal1(segPlaneX, segPlaneY, segPlaneZ);
  Vector oNormal2(nx, ny, nz);
  Vector oline = oNormal1 ^ oNormal2;
  real8 lenDeduc;
  Vector oPointToLeading;
  if (oline.Length() < dTol && (home->cycle % freq == 0)) {
    Vector odis(dx, dy, dz);
    odis.Normalize();
    lenDeduc = ovelo * odis;
    oPointToLeading = ovelo - ovelo * lenDeduc;
  } else if (home->cycle % freq != 0)
    return (SPLIT_FAILED);
  else {
    lenDeduc = oline * ovelo;
    ovelo = ovelo - oline * lenDeduc;
    lenDeduc = ovelo * oNormal1;
    oPointToLeading = ovelo - oNormal1 * lenDeduc;
  }
  oPointToLeading.Normalize();

  DissociateBurgersVector(home, bx, by, bz, &bx1, &by1, &bz1, &bx2, &by2, &bz2,
                          segPlaneX, segPlaneY, segPlaneZ, dx, dy, dz,
                          oPointToLeading);

  lastNode = nbrNode2;
  curNode = node;
  curNode2 = nbrNode2; // == poend2
  do {
    newNode = GetNewNativeNode(home);
    FreeNodeArms(newNode);
    newNode->constraint = UNCONSTRAINED;
    newNode->dNx = curNode->dNx;
    newNode->dNy = curNode->dNy;
    newNode->dNz = curNode->dNz;

    newNode->oldx = curNode->x;
    newNode->oldy = curNode->y;
    newNode->oldz = curNode->z;
    newNode->csState = ON_CROSS_SLIP_PLANE;

    pos[X] = curNode->x + lenDissociate * oPointToLeading.GetX();
    pos[Y] = curNode->y + lenDissociate * oPointToLeading.GetY();
    pos[Z] = curNode->z + lenDissociate * oPointToLeading.GetZ();
    // set the new node position
    newNode->x = pos[X];
    newNode->y = pos[Y];
    newNode->z = pos[Z];
    if (globalOp) {
      AddOp(home, SPLIT_NODE, node->myTag.domainID, node->myTag.index,
            newNode->myTag.domainID, newNode->myTag.index, -1,
            -1,            /* 3rd node tag unneeded */
            0.0, 0.0, 0.0, /* burgers vector unneeded*/
            newNode->x, newNode->y, newNode->z, 0.0, 0.0,
            0.0); // normal not needed
    }
    if (lastNode == poend2 && poend2->constraint >=SURFACE_NODE) {
      poend2->csState = ON_CROSS_SLIP_PLANE;
      Node_t *newNode1 = GetNewNativeNode(home);
      FreeNodeArms(newNode1);
      newNode1->constraint = SURFACE_NODE;
      newNode1->dNx = poend2->dNx;
      newNode1->dNy = poend2->dNy;
      newNode1->dNz = poend2->dNz;
      Vector oSur(poend2->dNx, poend2->dNy, poend2->dNz);
      lenDeduc = oSur * oPointToLeading;
      oPointToLeading = oPointToLeading - oSur * lenDeduc;
      oPointToLeading.Normalize();
      newNode1->oldx = poend2->x;
      newNode1->oldy = poend2->y;
      newNode1->oldz = poend2->z;
      newNode1->csState = ON_CROSS_SLIP_PLANE;
      newNode1->x = poend2->x + lenDissociate * oPointToLeading.GetX();
      newNode1->y = poend2->y + lenDissociate * oPointToLeading.GetY();
      newNode1->z = poend2->z + lenDissociate * oPointToLeading.GetZ();
      if (globalOp) {
        AddOp(home, SPLIT_NODE, node->myTag.domainID, node->myTag.index,
              newNode1->myTag.domainID, newNode1->myTag.index, -1,
              -1,            /* 3rd node tag unneeded */
              0.0, 0.0, 0.0, /* burgers vector unneeded*/
              newNode1->x, newNode1->y, newNode1->z, 0.0, 0.0,
              0.0); // normal not needed
      }
      InsertArm(home, newNode, &newNode1->myTag, bx1, by1, bz1, segPlaneX,
                segPlaneY, segPlaneZ, iChainID, globalOp);
      InsertArm(home, newNode1, &newNode->myTag, -bx1, -by1, -bz1, segPlaneX,
                segPlaneY, segPlaneZ, iChainID, globalOp);
    } else {
      InsertArm(home, newNode, &lastNode->myTag, bx1, by1, bz1, segPlaneX,
                segPlaneY, segPlaneZ, iChainID, globalOp);
      InsertArm(home, lastNode, &newNode->myTag, -bx1, -by1, -bz1, segPlaneX,
                segPlaneY, segPlaneZ, iChainID, globalOp);
    }

    nbrNode1 = GetNeighborNode(home, curNode, 0);
    nbrNode2 = GetNeighborNode(home, curNode, 1);
    if (nbrNode2 != curNode2) {
      nbrNode1 = nbrNode2;
      nbrNode2 = curNode2;
    }
    curNode2 = curNode;
    curNode = nbrNode1;
    lastNode = newNode;
  } while (curNode != poend1);

  if (poend1->constraint >=SURFACE_NODE) {
    poend1->csState = ON_CROSS_SLIP_PLANE;
    Node_t *newNode1 = GetNewNativeNode(home);
    FreeNodeArms(newNode1);
    newNode1->constraint = SURFACE_NODE;
    newNode1->dNx = poend1->dNx;
    newNode1->dNy = poend1->dNy;
    newNode1->dNz = poend1->dNz;
    Vector oSur(poend1->dNx, poend1->dNy, poend1->dNz);
    lenDeduc = oSur * oPointToLeading;
    oPointToLeading = oPointToLeading - oSur * lenDeduc;
    oPointToLeading.Normalize();
    newNode1->oldx = poend1->x;
    newNode1->oldy = poend1->y;
    newNode1->oldz = poend1->z;
    newNode1->csState = ON_CROSS_SLIP_PLANE;
    newNode1->x = poend1->x + lenDissociate * oPointToLeading.GetX();
    newNode1->y = poend1->y + lenDissociate * oPointToLeading.GetY();
    newNode1->z = poend1->z + lenDissociate * oPointToLeading.GetZ();
    if (globalOp) {
      AddOp(home, SPLIT_NODE, node->myTag.domainID, node->myTag.index,
            newNode1->myTag.domainID, newNode1->myTag.index, -1,
            -1,            /* 3rd node tag unneeded */
            0.0, 0.0, 0.0, /* burgers vector unneeded*/
            newNode1->x, newNode1->y, newNode1->z, 0.0, 0.0,
            0.0); // normal not needed
    }
    InsertArm(home, newNode, &newNode1->myTag, -bx1, -by1, -bz1, segPlaneX,
              segPlaneY, segPlaneZ, iChainID, globalOp);
    InsertArm(home, newNode1, &newNode->myTag, bx1, by1, bz1, segPlaneX,
              segPlaneY, segPlaneZ, iChainID, globalOp);
  } else {
    InsertArm(home, newNode, &poend1->myTag, -bx1, -by1, -bz1, segPlaneX,
              segPlaneY, segPlaneZ, iChainID, globalOp);
    InsertArm(home, poend1, &newNode->myTag, bx1, by1, bz1, segPlaneX,
              segPlaneY, segPlaneZ, iChainID, globalOp);
  }

  curNode = node;
  curNode2 = poend2; // == poend2
  do {
    curNode->csState = ON_CROSS_SLIP_PLANE;
    nbrNode1 = GetNeighborNode(home, curNode, 0);
    nbrNode2 = GetNeighborNode(home, curNode, 1);
    if (nbrNode2 != curNode2) {
      nbrNode1 = nbrNode2;
      nbrNode2 = curNode2;
    }
    ChangeArmBurg(home, curNode, &nbrNode2->myTag, bx2, by2, bz2, segPlaneX,
                  segPlaneY, segPlaneZ, globalOp, DEL_SEG_NONE, false);
    ChangeArmBurg(home, nbrNode2, &curNode->myTag, -bx2, -by2, -bz2, segPlaneX,
                  segPlaneY, segPlaneZ, globalOp, DEL_SEG_NONE, false);
    curNode2 = curNode;
    curNode = nbrNode1;
  } while (curNode != poend1);
  curNode = curNode2;
  ChangeArmBurg(home, curNode, &nbrNode1->myTag, -bx2, -by2, -bz2, segPlaneX,
                segPlaneY, segPlaneZ, globalOp, DEL_SEG_NONE, false);
  ChangeArmBurg(home, nbrNode1, &curNode->myTag, bx2, by2, bz2, segPlaneX,
                segPlaneY, segPlaneZ, globalOp, DEL_SEG_NONE, false);
}

/* This function dissociate a burgers vector according to the Thompson
 * tetrahedron Input bx, return the dissociation result(two burger vectors)
 */
void DissociateBurgersVector(Home_t *home, const real8 bx, const real8 by,
                             const real8 bz, real8 *bx1, real8 *by1, real8 *bz1,
                             real8 *bx2, real8 *by2, real8 *bz2, const real8 nx,
                             const real8 ny, const real8 nz, const real8 dx,
                             const real8 dy, const real8 dz,
                             const Vector oPointToLeading) {
  real8 dTolerance = 1.0E-1;
  Vector oBurgersVector(bx, by, bz);
  Vector oNormal(nx, ny, nz);
  Vector oline(dx, dy, dz);
  oline.Normalize();
  Vector oSecond = oNormal ^ oline;
  Vector oReferenceBurgers1(0.0, 0.0, 0.0);
  Vector oReferenceBurgers2(0.0, 0.0, 0.0);
  Vector oReferenceBurgers3(0.0, 0.0, 0.0);
  // slip plane 1,1,1
  oReferenceBurgers1.Set(1.0, -1.0, 0.0);
  oReferenceBurgers1.Normalize();
  oReferenceBurgers2.Set(0.0, 1.0, -1.0);
  oReferenceBurgers2.Normalize();
  oReferenceBurgers3.Set(-1.0, 0.0, 1.0);
  oReferenceBurgers3.Normalize();
  Vector oPartial(2.0 / 3.0 / sqrt(2.0), 1.0 / 3.0 / sqrt(2.0),
                  1.0 / 3.0 / sqrt(2.0));

  if (oReferenceBurgers1.IsSame(oBurgersVector, dTolerance)) {
    *bx1 = oPartial.GetX();
    *by1 = -oPartial.GetY();
    *bz1 = -oPartial.GetZ();
    *bx2 = oPartial.GetY();
    *by2 = -oPartial.GetX();
    *bz2 = oPartial.GetZ();
  } else if (oReferenceBurgers1.IsOpposite(oBurgersVector, dTolerance)) {
    *bx1 = -oPartial.GetY();
    *by1 = oPartial.GetX();
    *bz1 = -oPartial.GetZ();
    *bx2 = -oPartial.GetX();
    *by2 = oPartial.GetY();
    *bz2 = oPartial.GetZ();
  } else if (oReferenceBurgers2.IsSame(oBurgersVector, dTolerance)) {
    *bx1 = -oPartial.GetY();
    *by1 = oPartial.GetX();
    *bz1 = -oPartial.GetZ();
    *bx2 = oPartial.GetZ();
    *by2 = oPartial.GetY();
    *bz2 = -oPartial.GetX();
  } else if (oReferenceBurgers2.IsOpposite(oBurgersVector, dTolerance)) {
    *bx1 = -oPartial.GetZ();
    *by1 = -oPartial.GetY();
    *bz1 = oPartial.GetX();
    *bx2 = oPartial.GetY();
    *by2 = -oPartial.GetX();
    *bz2 = oPartial.GetZ();
  } else if (oReferenceBurgers3.IsSame(oBurgersVector, dTolerance)) {
    *bx1 = -oPartial.GetY();
    *by1 = -oPartial.GetZ();
    *bz1 = oPartial.GetX();
    *bx2 = -oPartial.GetX();
    *by2 = oPartial.GetY();
    *bz2 = oPartial.GetZ();
  } else if (oReferenceBurgers3.IsOpposite(oBurgersVector, dTolerance)) {
    *bx1 = oPartial.GetX();
    *by1 = -oPartial.GetY();
    *bz1 = -oPartial.GetZ();
    *bx2 = oPartial.GetZ();
    *by2 = oPartial.GetY();
    *bz2 = -oPartial.GetX();
  }
  if (oPointToLeading * oSecond >= 1 - dTolerance) {
    real8 tmp;
    tmp = *bx1;
    *bx1 = *bx2;
    *bx2 = tmp;
    tmp = *by1;
    *by1 = *by2;
    *by2 = tmp;
    tmp = *bz1;
    *bz1 = *bz2;
    *bz2 = tmp;
  } else if (oPointToLeading * oSecond <= -1 + dTolerance)
    ;
}
/*
  real8 dTolerance = 1.0E-1;
  Vector oBurgersVector(bx, by, bz);
  Vector oNormal(nx, ny, nz);
  Vector oline(dx, dy, dz);
  oline.Normalize();
  Vector oSecond = oNormal ^ oline;
  Vector oReferenceBurgers1(0.0, 0.0, 0.0);
  Vector oReferenceBurgers2(0.0, 0.0, 0.0);
  Vector oReferenceBurgers3(0.0, 0.0, 0.0);
  // slip plane -1,1,-1
  oReferenceBurgers1.Set(-1.0, -1.0, 0.0);
  oReferenceBurgers1.Normalize();
  oReferenceBurgers2.Set(0.0, 1.0, 1.0);
  oReferenceBurgers2.Normalize();
  oReferenceBurgers3.Set(-1.0, 0.0, 1.0);
  oReferenceBurgers3.Normalize();
  Vector oPartial(2.0 / 3.0 / sqrt(2.0), 1.0 / 3.0 / sqrt(2.0),
                  1.0 / 3.0 / sqrt(2.0));

  if (oReferenceBurgers1.IsSame(oBurgersVector, dTolerance)) {
    *bx1 = -oPartial.GetX();
    *by1 = -oPartial.GetY();
    *bz1 = oPartial.GetZ();
    *bx2 = -oPartial.GetY();
    *by2 = -oPartial.GetX();
    *bz2 = -oPartial.GetZ();
  } else if (oReferenceBurgers1.IsOpposite(oBurgersVector, dTolerance)) {
    *bx1 = oPartial.GetY();
    *by1 = oPartial.GetX();
    *bz1 = oPartial.GetZ();
    *bx2 = oPartial.GetX();
    *by2 = oPartial.GetY();
    *bz2 = -oPartial.GetZ();
  } else if (oReferenceBurgers2.IsSame(oBurgersVector, dTolerance)) {
    *bx1 = oPartial.GetY();
    *by1 = oPartial.GetX();
    *bz1 = oPartial.GetZ();
    *bx2 = -oPartial.GetZ();
    *by2 = oPartial.GetY();
    *bz2 = oPartial.GetX();
  } else if (oReferenceBurgers2.IsOpposite(oBurgersVector, dTolerance)) {
    *bx1 = oPartial.GetZ();
    *by1 = -oPartial.GetY();
    *bz1 = -oPartial.GetX();
    *bx2 = -oPartial.GetY();
    *by2 = -oPartial.GetX();
    *bz2 = -oPartial.GetZ();
  } else if (oReferenceBurgers3.IsSame(oBurgersVector, dTolerance)) {
    *bx1 = -oPartial.GetX();
    *by1 = -oPartial.GetZ();
    *bz1 = oPartial.GetY();
    *bx2 = -oPartial.GetZ();
    *by2 = oPartial.GetY();
    *bz2 = oPartial.GetX();
  } else if (oReferenceBurgers3.IsOpposite(oBurgersVector, dTolerance)) {
    *bx1 = oPartial.GetZ();
    *by1 = -oPartial.GetY();
    *bz1 = -oPartial.GetX();
    *bx2 = oPartial.GetX();
    *by2 = oPartial.GetY();
    *bz2 = -oPartial.GetZ();
  }
  if (oPointToLeading * oSecond >= 1 - dTolerance) {
    real8 tmp;
    tmp = *bx1;
    *bx1 = *bx2;
    *bx2 = tmp;
    tmp = *by1;
    *by1 = *by2;
    *by2 = tmp;
    tmp = *bz1;
    *bz1 = *bz2;
    *bz2 = tmp;
  } else if (oPointToLeading * oSecond <= -1 + dTolerance)
    ;
}*/
