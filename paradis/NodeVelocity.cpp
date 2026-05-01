/**************************************************************************
 *
 *      Module:       NodeVelocity.c
 *      Description:  Contains functions to control setting nodal
 *                    velocities, generating velocity statistics,
 *                    and setting the timestep.
 *
 ***************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>

#include "Home.h"
#include "Mobility.h"
#include "Util.h"

#ifdef PARALLEL
#include "mpi.h"
#endif

/*
 *     Define the types of velocity statistics to be accumulated.
 *     NOTE: V_MAXSTATS should always be the final item in the
 *     enumerated list below.
 */
typedef enum {
  V_NODE_COUNT = 0,
  V_AVERAGE_X,
  V_AVERAGE_Y,
  V_AVERAGE_Z,
  V_VAR,
  V_MAXSTATS
} VStatTypes_t;

/*-------------------------------------------------------------------------
 *
 *      Function:    GetVelocityStatistics
 *      Description: If gathering of velocity statistics is enabled,
 *                   gather the statistics (defined by the VStatTypes_t
 *                   above)
 *
 *------------------------------------------------------------------------*/
void GetVelocityStatistics(Home_t *home) {
#ifdef VEL_STATISTICS
  int i, nodeCount;
  real8 velStatsLocal[V_MAXSTATS], velStatsGlobal[V_MAXSTATS];
  real8 v2, vx, vy, vz;
  real8 v2sq, vAveragesq, vStDev;
  real8 vAverage, vAveragex, vAveragey, vAveragez;
  Param_t *param;
  Node_t *node;

  param = home->param;

  /*
   *      Zero out some stuff before starting
   */
  for (i = 0; i < V_MAXSTATS; i++) {
    velStatsLocal[i] = 0.0;
    velStatsGlobal[i] = 0.0;
  }

  /*
   *      Loop over all the nodes, accumulating all the necessary data
   */
  for (i = 0; i < home->newNodeKeyPtr; i++) {

    if ((node = home->nodeKeys[i]) == (Node_t *)NULL) {
      continue;
    }

    vx = node->vX;
    vy = node->vY;
    vz = node->vZ;

    v2 = vx * vx + vy * vy + vz * vz;

    /*
     *        If we're gathering statistics, accumulate the necessary data.
     *        Otherwise just find the highest nodal velocity.
     */
    velStatsLocal[V_AVERAGE_X] += vx;
    velStatsLocal[V_AVERAGE_Y] += vy;
    velStatsLocal[V_AVERAGE_Z] += vz;
    velStatsLocal[V_VAR] += v2;
    velStatsLocal[V_NODE_COUNT]++;
  }

  if (velStatsLocal[V_NODE_COUNT] > 0) {

    MPI_Allreduce(velStatsLocal, velStatsGlobal, V_MAXSTATS, MPI_DOUBLE,
                  MPI_SUM, MPI_COMM_WORLD);

    nodeCount = velStatsGlobal[V_NODE_COUNT];

    vAveragex = velStatsGlobal[V_AVERAGE_X] / nodeCount;
    vAveragey = velStatsGlobal[V_AVERAGE_Y] / nodeCount;
    vAveragez = velStatsGlobal[V_AVERAGE_Z] / nodeCount;

    vAveragesq =
        vAveragex * vAveragex + vAveragey * vAveragey + vAveragez * vAveragez;

    vStDev = sqrt(velStatsGlobal[V_VAR] / nodeCount) - vAveragesq;

    param->vStDev = vStDev;
    param->vAverage = sqrt(vAveragesq);
  }

#endif /* ifdef VEL_STATISTICS */

  return;
}

void AddFrictionForces(Home_t *poHome) {
  // this function is disabled for now
  return;
  // if there is no friction in the system, return
  double dTolerance = 1.0E-3;
  if ((poHome->param->FrictionStress <= dTolerance) &&
      (poHome->param->SolidSolutionStrength <= dTolerance) &&
      (poHome->param->KWStress <= dTolerance)) {
    return;
  }
  unsigned int i = 0;
  unsigned int j = 0;
  Node_t *poNode = NULL;
  Node_t *poNeighbour = NULL;
  Vector oLineDirection;
  Vector oBurgersVector;
  double dLength = 0.0;
  double dFrictionStress = 0.0;
  Vector oFrictionForce;
  Vector oArmForce;
  Vector oArmForceDirection;
  double dFrictionForceMagnitude = 0.0;
  Polyhedron *poPrecipitate = NULL;
  // loop over all the arms and subtraction the friction forces from the total
  // arm force
  for (i = 0; i < poHome->newNodeKeyPtr; i++) {
    poNode = poHome->nodeKeys[i];
    if (poNode == NULL)
      continue;
    dFrictionStress = 0.0;
    if (poHome->param->FrictionStress > 0.0)
      dFrictionStress += poHome->param->FrictionStress;
    if (poHome->param->KWStress > 0.0) {
      if (poHome->poPrecipitateServer->IsNodeInsideAPrecipitate(
              poNode, poPrecipitate)) {
        dFrictionStress += poHome->param->KWStress;
      } else if (poHome->param->SolidSolutionStrength > 0.0) {
        dFrictionStress += poHome->param->SolidSolutionStrength;
      }
    } else if (poHome->param->SolidSolutionStrength > 0.0) {
      dFrictionStress += poHome->param->SolidSolutionStrength;
    }

    for (j = 0; j < poNode->numNbrs; j++) {
      poNeighbour = GetNeighborNode(poHome, poNode, j);
      oLineDirection.SetByPoints(
          Point(poNode->x, poNode->y, poNode->z),
          Point(poNeighbour->x, poNeighbour->y, poNeighbour->z));
      dLength = oLineDirection.Length();
      oLineDirection.Normalize();
      oBurgersVector.Set(poNode->burgX[j], poNode->burgY[j], poNode->burgZ[j]);
      // The friction force is given by f_i = -tau b_j l s_i s_j
      // where tau is the friction stress, b_j are the Burgers vector
      // components, l is the segment length and s_i are the components of the
      // segment displacement normal to the dislocation segment. The negative
      // sign means that the friction opposes the motion. To get the friction
      // force, project the arm force on a direction normal to the dislocation
      // segment, normalize the resulting vector, this is the normal
      // displacement vector s_i. Given the segment length, friction stress and
      // the segment's Burges vector, get the friction force and subtract it
      // from the total force. If the friction force is greater than the total
      // force, zero out the total force
      oArmForce.Set(poNode->armfx[j], poNode->armfy[j],
                    poNode->armfz[j]); // junjie: should multiply by 2, because
                                       // armf is half of the segment force
      oArmForce = oArmForce - oLineDirection * (oArmForce * oLineDirection);
      oArmForceDirection = oArmForce;
      oArmForceDirection.Normalize();
      dFrictionForceMagnitude = fabs(dLength * dFrictionStress *
                                     (oBurgersVector * oArmForceDirection));
      if (dFrictionForceMagnitude > oArmForce.Length()) {
        poNode->armfx[j] = 0.0;
        poNode->armfy[j] = 0.0;
        poNode->armfz[j] = 0.0;
        continue;
      }
      oFrictionForce = oArmForceDirection * dFrictionForceMagnitude;
      poNode->armfx[j] -= oFrictionForce.GetX();
      poNode->armfy[j] -= oFrictionForce.GetY();
      poNode->armfz[j] -= oFrictionForce.GetZ();
    }
  }
  // update the nodal forces
  for (i = 0; i < poHome->newNodeKeyPtr; i++) {
    poNode = poHome->nodeKeys[i];
    if (poNode == NULL)
      continue;
    poNode->fX = 0.0;
    poNode->fY = 0.0;
    poNode->fZ = 0.0;
    for (j = 0; j < poNode->numNbrs; j++) {
      poNode->fX += poNode->armfx[j];
      poNode->fY += poNode->armfy[j];
      poNode->fZ += poNode->armfz[j];
    }
  }
}

// junjie
void AddPeierlsStress(Home_t *poHome) {
  // if the Peierls Stress is disabled, return
  double dAngle = PI / 180.0 * 10.0; // tolerance for screw dislocations
  if (poHome->param->EnablePeierls <= 0) {
    return;
  }
  unsigned int i = 0;
  unsigned int j = 0;
  Node_t *poNode = NULL;
  Node_t *poNeighbour = NULL;
  Vector oLineDirection;
  Vector oBurgersVector;
  double dLength = 0.0;
  Vector oFrictionForce;
  Vector oArmForce;
  Vector oArmForceDirection;
  Vector oNormal;
  Vector oBurgersVectorDirection;
  double dFrictionForceMagnitude = 0.0;
  double dFrictionStress = 0.0;
  // loop over all the arms and subtract the Peierls forces from the total arm
  // force
  for (i = 0; i < poHome->newNodeKeyPtr; i++) {
    poNode = poHome->nodeKeys[i];
    if (poNode == NULL)
      continue;
    for (j = 0; j < poNode->numNbrs; j++) {
      poNeighbour = GetNeighborNode(poHome, poNode, j);
      oLineDirection.SetByPoints(
          Point(poNode->x, poNode->y, poNode->z),
          Point(poNeighbour->x, poNeighbour->y, poNeighbour->z));
      dLength = oLineDirection.Length();
      oLineDirection.Normalize();
      oBurgersVector.Set(poNode->burgX[j], poNode->burgY[j], poNode->burgZ[j]);
      oBurgersVectorDirection = oBurgersVector;
      oBurgersVectorDirection.Normalize();
      oNormal.Set(poNode->nx[j], poNode->ny[j], poNode->nz[j]);
      // The Peierls force is given by f_i = -tau l b s_i
      // where tau is the Peierls stress for this component, l is
      // the segment length and s_i are the components of the segment
      // displacement normal to the dislocation segment. The negative sign means
      // that the Peierls Stress opposes the motion. To get the Peierls force,
      // project the arm force on to the slip plane of this segment, and then,
      // project the projection to the direction normal to the dislocation
      // segment, normalize the resulting vector, this is the normal
      // displacement vector s_i. Given the segment length, segment
      // component(screw or others), get the Peierls force and subtract it from
      // the total force. If the Peierls force is greater than the total force,
      // zero out the total force
      oArmForce.Set(2 * poNode->armfx[j], 2 * poNode->armfy[j],
                    2 * poNode->armfz[j]);
      oArmForce = oArmForce - oNormal * (oArmForce * oNormal);
      oArmForce = oArmForce - oLineDirection * (oArmForce * oLineDirection);
      oArmForceDirection = oArmForce;
      oArmForceDirection.Normalize();
      // check the component of this segment(screw or not screw)
      if (fabs(oLineDirection * oBurgersVectorDirection) > cos(dAngle)) {
        dFrictionStress = poHome->param->PeierlsScrew;
      } else {
        dFrictionStress = poHome->param->PeierlsOthers;
      }
      dFrictionForceMagnitude =
          fabs(dLength * dFrictionStress * (oBurgersVectorDirection.Length()));
      if (dFrictionForceMagnitude > oArmForce.Length()) {
        poNode->armfx[j] = 0.0;
        poNode->armfy[j] = 0.0;
        poNode->armfz[j] = 0.0;
        continue;
      }
      oFrictionForce = oArmForceDirection * (dFrictionForceMagnitude / 2);
      poNode->armfx[j] -= oFrictionForce.GetX();
      poNode->armfy[j] -= oFrictionForce.GetY();
      poNode->armfz[j] -= oFrictionForce.GetZ();
    }
  }
  // update the nodal forces
  for (i = 0; i < poHome->newNodeKeyPtr; i++) {
    poNode = poHome->nodeKeys[i];
    if (poNode == NULL)
      continue;
    poNode->fX = 0.0;
    poNode->fY = 0.0;
    poNode->fZ = 0.0;
    for (j = 0; j < poNode->numNbrs; j++) {
      poNode->fX += poNode->armfx[j];
      poNode->fY += poNode->armfy[j];
      poNode->fZ += poNode->armfz[j];
    }
  }
}
// junjie

/*-------------------------------------------------------------------------
 *
 *      Function:    CalcNodeVelocities
 *      Description: Driver function that will invoke the appropriate
 *                   mobility function to update the velocity of every
 *                   native node, then apply the velocity cutoff if
 *                   applicable.
 *
 *      Arguments:
 *          doAll      Flag indicating if ALL nodes are to have
 *                     veolcity recalculated, or just those nodes
 *                     that have had their forces updated.
 *
 *      Returns:  0 if velocity was successfully calculated for all
 *                  nodes
 *                1 if the mobility functions were unable to converge
 *                  on a velocity for one or more nodes.
 *
 *------------------------------------------------------------------------*/
int CalcNodeVelocities(Home_t *home) {
  // Add the friction forces to the arm and node forces
  AddFrictionForces(home);
  // junjie: add the Peierls stress to the arm and node forces
  AddPeierlsStress(home);
  // junjie

  int domainMobError;
  Node_t *node;
  Param_t *param;

  param = home->param;
  domainMobError = 0;

  int i, nodeMobError;
  Node_t *poNode;

  for (i = 0; i < home->newNodeKeyPtr; i++) {
    /*
     *              If we encountered a mobility error on a previous node,
     *              we'll probably be cutting the timestep, so don't
     *              bother restting the velocity of any subsequent nodes.
     */
    if (domainMobError != 0) {
      break;
    }

    if ((poNode = home->nodeKeys[i]) == NULL) {
      continue;
    }

    /*
     *              We set a pointer to the appropriate mobility function
     *              during initialization, so just invoke the function now.
     *              Need two error flags.  One to indicate if the current node
     *              had a mobility error, the other (which is returned to the
     *              caller) indicates if there were *any* nodes with mobility
     *              errors in the domain.
     */
    nodeMobError = param->mobilityFunc(home, poNode);
    domainMobError |= nodeMobError;
  }

  GetVelocityStatistics(home);

#if PARALLEL
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  return (domainMobError);
}
