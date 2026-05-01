/**************************************************************************
 *
 *  Function    : Mobility_FCC_0
 *  Author      : Wei Cai, Seok-Woo Lee (updated 07/14/09)
 *  Description : Generic Mobility Law of FCC metals
 *                Each line has a glide plane from input
 *                and it never changes
 *                If the plane normal is not of {111} type, dislocation
 *                motion is constrained along line direction
 *                If node flag == 7, node velocity is zero
 *
 *  Returns:  0 on success
 *            1 if velcoity could not be determined
 *
 ***************************************************************************/

#include "Home.h"
#include "Util.h"
#include "Mobility.h"
#include <stdio.h>
#include <math.h>
#include "ParadisSurface.h"
#include "ParadisCrossSlipServer.h"

#ifdef PARALLEL
#include "mpi.h"
#endif

int Mobility_FCC_0_Old(Home_t *home, Node_t *node) {
  int numNonZeroLenSegs = 0;
  Param_t *param;
  real8 VelxNode, VelyNode, VelzNode, Veldotlcr;
  int i, j, nc, nconstraint, nlc;
  int n_x, n_y, n_z, s_timestep;
  real8 normX[100], normY[100], normZ[100], normx[100], normy[100], normz[100];
  real8 lineX[100], lineY[100], lineZ[100];
  real8 a, b;
  real8 dx, dy, dz, lx, ly, lz, lr, LtimesB;
  real8 lcx, lcy, lcz, normdotlc;
  Node_t *nbr;
  real8 MobScrew, MobEdge, Mob;
  real8 bx, by, bz, br, dangle;
  real8 nForce[3];
  real8 pos_mid[3], pos_stress0[3], pos_stress1[3];
  real8 dTemperature, xd, yd, zd;
  real8 c000, c001, c010, c100, c110, c101, c011, c111;
  param = home->param;

  MobScrew = param->MobScrew;
  MobEdge = param->MobEdge;

  nc = node->numNbrs;

  /*
   *  If node is 'pinned' in place by constraints, or the node has any arms
   *  with a burgers vector that has explicitly been set to be sessile (via
   *  control file inputs), the node may not be moved so just zero the velocity
   *  and return
   */
  if ((node->constraint == PINNED_NODE) || NodeHasSessileBurg(home, node)) {
    node->vX = 0.0;
    node->vY = 0.0;
    node->vZ = 0.0;
    return (0);
  }

  /*
   *  It's possible this function was called for a node which had only zero-
   *  length segments (during SplitSurfaceNodes() for example).  If that is
   *  the case, just set the velocity to zero and return.
   */
  for (i = 0; i < nc; i++) {
    if ((nbr = GetNeighborNode(home, node, i)) == NULL) {
      continue;
    }
    dx = node->x - nbr->x;
    dy = node->y - nbr->y;
    dz = node->z - nbr->z;
    if ((dx * dx + dy * dy + dz * dz) > 1.0e-12) {
      numNonZeroLenSegs++;
    }
  }

  if (numNonZeroLenSegs == 0) {
    node->vX = 0.0;
    node->vY = 0.0;
    node->vZ = 0.0;
    return (0);
  }

  /* copy glide plane constraints and determine line constraints */
  for (i = 0; i < nc; i++) {
    normX[i] = normx[i] = node->nx[i];
    normY[i] = normy[i] = node->ny[i];
    normZ[i] = normz[i] = node->nz[i];

    if ((fabs(fabs(node->nx[i]) - fabs(node->ny[i])) > FFACTOR_NORMAL) ||
        (fabs(fabs(node->ny[i]) - fabs(node->nz[i])) > FFACTOR_NORMAL)) {
      /* not {111} plane */
      if ((nbr = GetNeighborNode(home, node, i)) == (Node_t *)NULL) {
        Fatal("Neighbor not found at %s line %d\n", __FILE__, __LINE__);
      }
      lineX[i] = nbr->x - node->x;
      lineY[i] = nbr->y - node->y;
      lineZ[i] = nbr->z - node->z;
    } else { /* no line constraint */
      lineX[i] = lineY[i] = lineZ[i] = 0;
    }
  }

  /* normalize glide plane normal vectors and lc line vectors*/
  for (i = 0; i < nc; i++) {
    a = sqrt(normX[i] * normX[i] + normY[i] * normY[i] + normZ[i] * normZ[i]);
    b = sqrt(lineX[i] * lineX[i] + lineY[i] * lineY[i] + lineZ[i] * lineZ[i]);

    if (a > 0) {
      normX[i] /= a;
      normY[i] /= a;
      normZ[i] /= a;

      normx[i] /= a;
      normy[i] /= a;
      normz[i] /= a;
    }
    if (b > 0) {
      lineX[i] /= b;
      lineY[i] /= b;
      lineZ[i] /= b;
    }
  }

  /* Find independent glide constraints */
  nconstraint = nc;
  for (i = 0; i < nc; i++) {
    for (j = 0; j < i; j++) {
      Orthogonalize(normX + i, normY + i, normZ + i, normX[j], normY[j],
                    normZ[j]);
    }
    if ((normX[i] * normX[i] + normY[i] * normY[i] + normZ[i] * normZ[i]) <
        FFACTOR_ORTH) {
      normX[i] = normY[i] = normZ[i] = 0.0;
      nconstraint--;
    }
  }

  /* Find independent line constraints */
  nlc = 0;
  for (i = 0; i < nc; i++) {
    for (j = 0; j < i; j++) {
      Orthogonalize(lineX + i, lineY + i, lineZ + i, lineX[j], lineY[j],
                    lineZ[j]);
    }
    if ((lineX[i] * lineX[i] + lineY[i] * lineY[i] + lineZ[i] * lineZ[i]) <
        FFACTOR_ORTH) {
      lineX[i] = lineY[i] = lineZ[i] = 0.0;
    } else {
      nlc++;
    }
  }

  /* find total dislocation length times drag coefficent (LtimesB)*/
  double dTolerance = 1.0E-6;
  LtimesB = 0;
  for (j = 0; j < nc; j++) {
    nbr = GetNeighborNode(home, node, j);
    if (nbr == NULL) {
      continue;
    }
    dx = nbr->x - node->x;
    dy = nbr->y - node->y;
    dz = nbr->z - node->z;

    lr = sqrt(dx * dx + dy * dy + dz * dz);
    if (lr < dTolerance) {
      /* zero arm segment can happen after node split
       * it is OK to have plane normal vector == 0
       * Skip (do nothing)
       */
      continue;
    }
    lx = dx / lr;
    ly = dy / lr;
    lz = dz / lr;

    bx = node->burgX[j];
    by = node->burgY[j];
    bz = node->burgZ[j];

    br = sqrt(bx * bx + by * by + bz * bz);
    bx = bx / br;
    by = by / br;
    bz = bz / br;

    dangle = fabs(bx * lx + by * ly + bz * lz);

    // Yejun
    if (param->Ttype == 1 && param->stress_timestep > 0) {
      pos_mid[0] = 0.5 * node->x + 0.5 * nbr->x;
      pos_mid[1] = 0.5 * node->y + 0.5 * nbr->y;
      pos_mid[2] = 0.5 * node->z + 0.5 * nbr->z;

      GetInLocalCoordinates(param, pos_mid[0], pos_mid[1], pos_mid[2]);
      if (fabs(pos_mid[0]) > 0.5 * param->stress_Dim[0] ||
          fabs(pos_mid[1]) > 0.5 * param->stress_Dim[1] ||
          fabs(pos_mid[2]) > 0.5 * param->stress_Dim[2]) {
        dTemperature = param->TempK;
      } else {
        n_x = ceil((pos_mid[0] + 0.5 * param->stress_Dim[0]) /
                   param->stress_Dim[0] * (double)(param->GP_x - 1));
        n_y = ceil((pos_mid[1] + 0.5 * param->stress_Dim[1]) /
                   param->stress_Dim[1] * (double)(param->GP_y - 1));
        n_z = ceil((pos_mid[2] + 0.5 * param->stress_Dim[2]) /
                   param->stress_Dim[2] * (double)(param->GP_z - 1));

        n_x = fmax(n_x, 1);
        n_y = fmax(n_y, 1);
        n_z = fmax(n_z, 1);

        s_timestep = param->stress_timestep - 1;
        pos_stress0[0] = home->stress[n_x - 1][n_y - 1][n_z - 1][s_timestep].x;
        pos_stress0[1] = home->stress[n_x - 1][n_y - 1][n_z - 1][s_timestep].y;
        pos_stress0[2] = home->stress[n_x - 1][n_y - 1][n_z - 1][s_timestep].z;
        pos_stress1[0] = home->stress[n_x][n_y][n_z][s_timestep].x;
        pos_stress1[1] = home->stress[n_x][n_y][n_z][s_timestep].y;
        pos_stress1[2] = home->stress[n_x][n_y][n_z][s_timestep].z;

        // trilinear interpolation
        xd = (pos_mid[0] - pos_stress0[0]) / (pos_stress1[0] - pos_stress0[0]);
        yd = (pos_mid[1] - pos_stress0[1]) / (pos_stress1[1] - pos_stress0[1]);
        zd = (pos_mid[2] - pos_stress0[2]) / (pos_stress1[2] - pos_stress0[2]);
        c000 = home->stress[n_x - 1][n_y - 1][n_z - 1][s_timestep].Temp;
        c001 = home->stress[n_x - 1][n_y - 1][n_z][s_timestep].Temp;
        c010 = home->stress[n_x - 1][n_y][n_z - 1][s_timestep].Temp;
        c011 = home->stress[n_x - 1][n_y][n_z][s_timestep].Temp;
        c100 = home->stress[n_x][n_y - 1][n_z - 1][s_timestep].Temp;
        c101 = home->stress[n_x][n_y - 1][n_z][s_timestep].Temp;
        c110 = home->stress[n_x][n_y][n_z - 1][s_timestep].Temp;
        c111 = home->stress[n_x][n_y][n_z][s_timestep].Temp;

        dTemperature = ((c000 * (1 - xd) + c100 * xd) * (1 - yd) +
                        (c010 * (1 - xd) + c110 * xd) * yd) *
                           (1 - zd) +
                       ((c001 * (1 - xd) + c101 * xd) * (1 - yd) +
                        (c011 * (1 - xd) + c111 * xd) * yd) *
                           zd;
      }
      Mob = 1.0 / (param->MobM * dTemperature + param->MobB0);
    } else {
      Mob = MobEdge + (MobScrew - MobEdge) * dangle;
    }
    LtimesB = LtimesB + (lr / Mob);
  }
  LtimesB = LtimesB / 2.0;

  nForce[0] = node->fX;
  nForce[1] = node->fY;
  nForce[2] = node->fZ;

  /* Velocity is simply proportional to total force per unit length */
  VelxNode = nForce[0] / LtimesB;
  VelyNode = nForce[1] / LtimesB;
  VelzNode = nForce[2] / LtimesB;

  /* Orthogonalize with glide plane constraints */
  for (i = 0; i < nc; i++) {
    if ((normX[i] != 0) || (normY[i] != 0) || (normZ[i] != 0)) {
      Orthogonalize(&VelxNode, &VelyNode, &VelzNode, normX[i], normY[i],
                    normZ[i]);
    }
  }

  /* Any dislocation with glide plane not {111} type can only move along its
   * length This rule includes LC junction which is on {100} plane
   */
  if (nlc == 1) {
    for (i = 0; i < nc; i++) {
      if ((fabs(lineX[i]) > dTolerance) || (fabs(lineY[i]) > dTolerance) ||
          (fabs(lineZ[i]) > dTolerance)) {
        lcx = lineX[i];
        lcy = lineY[i];
        lcz = lineZ[i];
        break;
      }
    }

    /* project velocity along line */
    Veldotlcr = VelxNode * lcx + VelyNode * lcy + VelzNode * lcz;
    VelxNode = Veldotlcr * lcx;
    VelyNode = Veldotlcr * lcy;
    VelzNode = Veldotlcr * lcz;

    if (nconstraint >= 1) {
      // a few plane constraints and one line constraint
      for (i = 0; i < nc; i++) {
        normdotlc = normx[i] * lcx + normy[i] * lcy + normz[i] * lcz;
        if (fabs(normdotlc) > FFACTOR_ORTH) {
          /* set velocity to zero if line is not on every plane */
          VelxNode = 0.0;
          VelyNode = 0.0;
          VelzNode = 0.0;
          break;
        }
      }
    }
  } else if (nlc >= 2) {
    /* Velocity is zero when # of independnet lc constratins is equal to or more
     * than 2 */
    VelxNode = 0.0;
    VelyNode = 0.0;
    VelzNode = 0.0;
  }

  node->vX = VelxNode;
  node->vY = VelyNode;
  node->vZ = VelzNode;

  if (node->constraint >=SURFACE_NODE) {
    unsigned int iConstraintType = 0;
    Vector oConstraintVector;
    ParadisSurface::GetNodeDynamicConstraint(node, iConstraintType,
                                             oConstraintVector);
    Vector oSurfaceNormal = Vector(node->dNx, node->dNy, node->dNz);
    if (iConstraintType == 1) {
      Vector oTargetVelocityDirection = oConstraintVector ^ oSurfaceNormal;
      oTargetVelocityDirection.Normalize();
      Vector oVelocity = Vector(node->vX, node->vY, node->vZ);
      double dTemp = oVelocity * oTargetVelocityDirection;
      node->vX = dTemp * oTargetVelocityDirection.GetX();
      node->vY = dTemp * oTargetVelocityDirection.GetY();
      node->vZ = dTemp * oTargetVelocityDirection.GetZ();
    } else if (iConstraintType == 2) {
      if (fabs(oSurfaceNormal * oConstraintVector) > dTolerance) {
        // velocity direction has a component in normal to the surface, which is
        // inconsistent with the surface motion constraint, the node should not
        // move at all
        node->vX = 0.0;
        node->vY = 0.0;
        node->vZ = 0.0;
      }
    } else {
      node->vX = 0.0;
      node->vY = 0.0;
      node->vZ = 0.0;
    }
  }

  return (0);
}

int Mobility_FCC_Glide(Home_t *poHome, Node_t *poNode) {
  // if node is pinned or has sessile bugers vector, set the velocity to zero
  // and return
  if (poNode->constraint == PINNED_NODE) {
    poNode->vX = 0.0;
    poNode->vY = 0.0;
    poNode->vZ = 0.0;
    return 0;
  }

  if (NodeHasSessileBurg(poHome, poNode)) {
    poNode->vX = 0.0;
    poNode->vY = 0.0;
    poNode->vZ = 0.0;
    return 0;
  }

  unsigned int iConstraintType = 0;
  Vector oConstraintVector;
  Vector armNormal;
  ParadisSurface::GetNodeDynamicConstraint(poNode, iConstraintType,
                                           oConstraintVector);
  // node is dynamically fixed, set the velocity to zero and return
  if (iConstraintType == 3) {
    poNode->vX = 0.0;
    poNode->vY = 0.0;
    poNode->vZ = 0.0;
    return 0;
  }

  // since this is FCC mobility, nodes that are not on {111} slip planes cannot
  // move, to be more accurate, they can only move along their arms, but this is
  // disabled for now
  if (iConstraintType == 1) {
    if (!ParadisCrossSlipServer::Is111Vector(oConstraintVector)) {
      poNode->vX = 0.0;
      poNode->vY = 0.0;
      poNode->vZ = 0.0;
      return 0;
    }
  }
  // junjie: since this is FCC mobility, intersection nodes that are not moving
  // along <110>
  // slip directions cannot move
  if (iConstraintType == 2) {
    if (!ParadisCrossSlipServer::Is110Vector(oConstraintVector)) {
      poNode->vX = 0.0;
      poNode->vY = 0.0;
      poNode->vZ = 0.0;
      return 0;
    }
    // junjie: nodes having a not {111} arm, cannot move
    for (int i = 0; i < poNode->numNbrs; ++i) {
      armNormal.Set(poNode->nx[i], poNode->ny[i], poNode->nz[i]);
      if (!ParadisCrossSlipServer::Is111Vector(armNormal)) {
        poNode->vX = 0.0;
        poNode->vY = 0.0;
        poNode->vZ = 0.0;
        return 0;
      }
    }
  }

  // store the old velocity
  Vector oOldVelocity(poNode->vX, poNode->vY, poNode->vZ);

  // calculate the node's mobility
  // find total dislocation length times drag coefficent (L x B)
  double dTolerance = 1.0E-6;
  double dLB = 0.0;
  unsigned int i = 0;
  Node_t *poNeighbour = NULL;
  double dX = 0.0;
  double dY = 0.0;
  double dZ = 0.0;
  double dLength = 0.0;
  double dBX = 0.0;
  double dBY = 0.0;
  double dBZ = 0.0;
  double dB = 0.0;
  double dProjection = 0.0;
  double dMobility = 0.0;
  double dScrewSegmentMobility = poHome->param->MobScrew;
  double dEdgeSegmentMobility = poHome->param->MobEdge;
  bool bIsLocked = true;
  int n_x, n_y, n_z, s_timestep;
  real8 pos_mid[3], pos_stress0[3], pos_stress1[3];
  real8 dTemperature, xd, yd, zd;
  real8 c000, c001, c010, c100, c110, c101, c011, c111;

  for (i = 0; i < poNode->numNbrs; i++) {
    poNeighbour = GetNeighborNode(poHome, poNode, i);
    if (poNeighbour == NULL) {
      continue;
    }
    if (poNode->constraint >=SURFACE_NODE) {
      if (poNeighbour->constraint >=SURFACE_NODE) {
        continue;
      }
    }

    dX = poNeighbour->x - poNode->x;
    dY = poNeighbour->y - poNode->y;
    dZ = poNeighbour->z - poNode->z;

    GetPrimaryImage(poHome->param, &dX, &dY, &dZ);

    // arm length
    dLength = sqrt(dX * dX + dY * dY + dZ * dZ);
    if (dLength < dTolerance) {
      // skip zero length segments
      continue;
    }

    // arm direction
    dX = dX / dLength;
    dY = dY / dLength;
    dZ = dZ / dLength;

    // burgers vector
    dBX = poNode->burgX[i];
    dBY = poNode->burgY[i];
    dBZ = poNode->burgZ[i];

    // burgers direction
    dB = sqrt(dBX * dBX + dBY * dBY + dBZ * dBZ);
    dBX = dBX / dB;
    dBY = dBY / dB;
    dBZ = dBZ / dB;
    dProjection = fabs(dBX * dX + dBY * dY + dBZ * dZ);

    // Yejun
    if (poHome->param->Ttype == 1 && poHome->param->stress_timestep > 0) {
      pos_mid[0] = 0.5 * poNode->x + 0.5 * poNeighbour->x;
      pos_mid[1] = 0.5 * poNode->y + 0.5 * poNeighbour->y;
      pos_mid[2] = 0.5 * poNode->z + 0.5 * poNeighbour->z;

      GetInLocalCoordinates(poHome->param, pos_mid[0], pos_mid[1], pos_mid[2]);
      if (fabs(pos_mid[0]) > 0.5 * poHome->param->stress_Dim[0] ||
          fabs(pos_mid[1]) > 0.5 * poHome->param->stress_Dim[1] ||
          fabs(pos_mid[2]) > 0.5 * poHome->param->stress_Dim[2]) {
        dTemperature = poHome->param->TempK;
      } else {
        n_x = ceil((pos_mid[0] + 0.5 * poHome->param->stress_Dim[0]) /
                   poHome->param->stress_Dim[0] *
                   (double)(poHome->param->GP_x - 1));
        n_y = ceil((pos_mid[1] + 0.5 * poHome->param->stress_Dim[1]) /
                   poHome->param->stress_Dim[1] *
                   (double)(poHome->param->GP_y - 1));
        n_z = ceil((pos_mid[2] + 0.5 * poHome->param->stress_Dim[2]) /
                   poHome->param->stress_Dim[2] *
                   (double)(poHome->param->GP_z - 1));

        n_x = fmax(n_x, 1);
        n_y = fmax(n_y, 1);
        n_z = fmax(n_z, 1);

        s_timestep = poHome->param->stress_timestep - 1;
        pos_stress0[0] =
            poHome->stress[n_x - 1][n_y - 1][n_z - 1][s_timestep].x;
        pos_stress0[1] =
            poHome->stress[n_x - 1][n_y - 1][n_z - 1][s_timestep].y;
        pos_stress0[2] =
            poHome->stress[n_x - 1][n_y - 1][n_z - 1][s_timestep].z;
        pos_stress1[0] = poHome->stress[n_x][n_y][n_z][s_timestep].x;
        pos_stress1[1] = poHome->stress[n_x][n_y][n_z][s_timestep].y;
        pos_stress1[2] = poHome->stress[n_x][n_y][n_z][s_timestep].z;

        // trilinear interpolation
        xd = (pos_mid[0] - pos_stress0[0]) / (pos_stress1[0] - pos_stress0[0]);
        yd = (pos_mid[1] - pos_stress0[1]) / (pos_stress1[1] - pos_stress0[1]);
        zd = (pos_mid[2] - pos_stress0[2]) / (pos_stress1[2] - pos_stress0[2]);
        c000 = poHome->stress[n_x - 1][n_y - 1][n_z - 1][s_timestep].Temp;
        c001 = poHome->stress[n_x - 1][n_y - 1][n_z][s_timestep].Temp;
        c010 = poHome->stress[n_x - 1][n_y][n_z - 1][s_timestep].Temp;
        c011 = poHome->stress[n_x - 1][n_y][n_z][s_timestep].Temp;
        c100 = poHome->stress[n_x][n_y - 1][n_z - 1][s_timestep].Temp;
        c101 = poHome->stress[n_x][n_y - 1][n_z][s_timestep].Temp;
        c110 = poHome->stress[n_x][n_y][n_z - 1][s_timestep].Temp;
        c111 = poHome->stress[n_x][n_y][n_z][s_timestep].Temp;

        dTemperature = ((c000 * (1 - xd) + c100 * xd) * (1 - yd) +
                        (c010 * (1 - xd) + c110 * xd) * yd) *
                           (1 - zd) +
                       ((c001 * (1 - xd) + c101 * xd) * (1 - yd) +
                        (c011 * (1 - xd) + c111 * xd) * yd) *
                           zd;
      }
      dMobility =
          1.0 / (poHome->param->MobM * dTemperature + poHome->param->MobB0);
    } else {
      dMobility = dEdgeSegmentMobility +
                  (dScrewSegmentMobility - dEdgeSegmentMobility) * dProjection;
    }
    dLB = dLB + (dLength / dMobility);
  }

  // if dLB is zero, zero the node's velocity
  if (dLB < dTolerance) {
    poNode->vX = 0.0;
    poNode->vY = 0.0;
    poNode->vZ = 0.0;
    return 0;
  }
  dLB = 0.5 * dLB;

  // velocity is proportional to total force per unit length
  Vector oVelocity = Vector(poNode->fX, poNode->fY, poNode->fZ);
  oVelocity = oVelocity * (1.0 / dLB);

  // now we have the velocity, project it according to the dynamic constraint,
  // in case of no dynamic constraints (constraint type 0), nothing needs to be
  // done
  if (iConstraintType == 1) {
    // node moves normal to the constraint vector, subtract the component of the
    // velocity vector in the direction of the constraint vector
    oVelocity = oVelocity - oConstraintVector * (oVelocity * oConstraintVector);
  } else if (iConstraintType == 2) {
    // node moves along the constraint vector, project the obtained velocity in
    // the direction of the constraint vector
    oVelocity = oConstraintVector * (oVelocity * oConstraintVector);
  }

  // if this node is a surface node, force it to move on the surface
  if (poNode->constraint >=SURFACE_NODE) {
    Vector oSurfaceNormal = Vector(poNode->dNx, poNode->dNy, poNode->dNz);
    if (iConstraintType == 1) {
      // node moves on a plane
      Vector oTargetVelocityDirection = oConstraintVector ^ oSurfaceNormal;
      oTargetVelocityDirection.Normalize();
      oVelocity =
          oTargetVelocityDirection * (oVelocity * oTargetVelocityDirection);
    } else if (iConstraintType == 2) {
      if (fabs(oSurfaceNormal * oConstraintVector) > dTolerance) {
        // velocity direction has a component in normal to the surface, which is
        // inconsistent with the surface motion constraint, the node should not
        // move at all
        oVelocity.Set(0.0, 0.0, 0.0);
      }
    } else {
      // node is free to move, just subtract the suface normal component of the
      // velocity vector from the velocity vector
      oVelocity = oVelocity - oSurfaceNormal * (oVelocity * oSurfaceNormal);
    }
  }

  // prevent node oscillations
  //     Vector oOldVelocityDirection = oOldVelocity.GetDirection();
  //     Vector oCurrentVelocityDirection = oVelocity.GetDirection();
  //     double dOppositeDirectionAngleCosine = 0.7;
  //     int MaximumVelocityDampingSteps = poHome->param->VelocityDampingSteps;
  //     double dVelocityDampingFactor = 1.0;
  //     if(poNode->VelocityDapmingSteps == 0)
  //     {
  //     	if(oOldVelocityDirection*oCurrentVelocityDirection <
  //     -dOppositeDirectionAngleCosine)
  //     	{
  //     		poNode->VelocityDapmingSteps =
  //     MaximumVelocityDampingSteps; 		dVelocityDampingFactor =
  //     pow(10.0,-poNode->VelocityDapmingSteps); 		oVelocity =
  //     oVelocity*dVelocityDampingFactor;
  //     	}
  //     }
  //     else
  //     {
  //     	poNode->VelocityDapmingSteps = poNode->VelocityDapmingSteps - 1;
  //     	dVelocityDampingFactor =
  //     pow(10.0,-poNode->VelocityDapmingSteps); 	oVelocity =
  //     oVelocity*dVelocityDampingFactor;
  //     }

  // finally, if the velocity vector is too high such that the time step will be
  // less than a preset hardcoded value, scale down the velocity vector
  //     double dVelocityMagnitude = oVelocity.Length();
  //     if(dVelocityMagnitude > 1.0)
  //     {
  //     	double dMinTimeStep = 5.0e-12;
  //     	double dTimeStep = poHome->param->rmax/dVelocityMagnitude;
  //     	if(dTimeStep < dMinTimeStep)
  //     	{
  //     		oVelocity = oVelocity*(dTimeStep/dMinTimeStep);
  //     	}
  //     }

  // set the node velocity
  poNode->vX = oVelocity.GetX();
  poNode->vY = oVelocity.GetY();
  poNode->vZ = oVelocity.GetZ();

  return 0;
}

#define MAX_SEG_PER_NODE 10

int Mobility_FCC_climb(Home_t *home, Node_t *node) {
  int i, j, numNbrs;
  int numNonZeroLenSegs = 0;
  real8 bDotb;
  real8 halfLen;
  real8 angle1, angle2;
  real8 costheta, costheta2, cosCritical1, cosCritical2;
  real8 Bline, Bscrew, Bedge, Bglide, Bclimb, Beclimb;
  real8 Bscrew2, Bedge2, Beclimb2;
  real8 invBscrew2, invBedge2;
  real8 BlmBecl;
  real8 eps = 1.0e-12;
  real8 segLen[MAX_SEG_PER_NODE];
  real8 burgList[MAX_SEG_PER_NODE][3];
  real8 segVector[MAX_SEG_PER_NODE][3];
  real8 nDir[MAX_SEG_PER_NODE][3];
  real8 lineDir[MAX_SEG_PER_NODE][3];
  real8 nForce[3], nVel[3];
  real8 Btotal[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
  real8 invBtotal[3][3];
  Param_t *param;

  param = home->param;

  Bscrew = 1.0 / param->MobScrew;
  Bedge = 1.0 / param->MobEdge;

  /*
   *      Climb is set this way to help with convergence and the fact that
   *      this is a glide restricted mobility
   */
  Beclimb = 10000.0 * Bedge;

  Bedge2 = Bedge * Bedge;
  Bscrew2 = Bscrew * Bscrew;
  Beclimb2 = Beclimb * Beclimb;

  Bline = 1.0 * MIN(Bscrew, Bedge);
  BlmBecl = Bline - Beclimb;

  invBscrew2 = 1.0 / (Bscrew * Bscrew);
  invBedge2 = 1.0 / (Bedge * Bedge);

  numNbrs = node->numNbrs;

  angle1 = 1.0;
  angle2 = 5.0;

  /*
   *      If node is 'pinned' in place by constraints, or the node has any arms
   *      with a burgers vector that has explicitly been set to be sessile (via
   *      control file inputs), the node may not be moved so just zero the
   *      velocity and return
   */
  if ((node->constraint == PINNED_NODE) || NodeHasSessileBurg(home, node)) {
    node->vX = 0.0;
    node->vY = 0.0;
    node->vZ = 0.0;
    return (0);
  }

  /*
   *      have to go through all node's segments and set up the
   *      appropriate array of glide planes.  While we're at it, might
   *      as well save some additional info we'll need later.
   */
  cosCritical1 = (M_PI * angle1 / 180.0);
  cosCritical2 = (M_PI * angle2 / 180.0);

  for (i = 0; i < numNbrs; i++) {
    real8 segLen2;
    Node_t *nbrNode;

    if ((nbrNode = GetNeighborNode(home, node, i)) == (Node_t *)NULL) {
      continue;
    }

    segVector[i][X] = node->x - nbrNode->x;
    segVector[i][Y] = node->y - nbrNode->y;
    segVector[i][Z] = node->z - nbrNode->z;

    /*
     *          Skip the zero length segments
     */
    segLen2 = DotProduct(segVector[i], segVector[i]);

    if ((segLen2) < eps) {
      segLen[i] = 0.0;
      continue;
    }

    numNonZeroLenSegs++;

    segLen[i] = sqrt(segLen2);

    /*
     *          Save some segment specific stuff for later rather
     */
    burgList[i][X] = node->burgX[i];
    burgList[i][Y] = node->burgY[i];
    burgList[i][Z] = node->burgZ[i];

    VECTOR_COPY(lineDir[i], segVector[i]);
    NormalizeVec(lineDir[i]);

    nDir[i][X] = node->nx[i];
    nDir[i][Y] = node->ny[i];
    nDir[i][Z] = node->nz[i];

    /*
     *          If this is one arm of a 2-node, see if it is a screw segment
     *
     *          The segment glide planes are calculated based on whether
     *          the segment is screw or not.  Two angles (angle1 and angle2)
     *          are defined which determine the 'screwness' of the segment.
     *          Any segment within <angle1> degrees of being screw is
     *          considered screw and its glide plane will simply be
     *          maintained.  A segment within <angle2> degrees of screw (but
     *          more than <angle1>) is considered to be in a 'mixed' region
     *          where a glide plane is calculated from a linear combination
     *          of the original screw glide plane and l cross b.  Glide planes
     *          for all other segments are simply calculated via l cross b.
     */
    if (numNbrs == 2) {

      bDotb = DotProduct(burgList[i], burgList[i]);
      costheta = DotProduct(lineDir[i], burgList[i]);
      costheta2 = (costheta * costheta) / bDotb;

      if (costheta2 <= (cosCritical1 * cosCritical1)) {
        /* not screw */
        NormalizedCrossVector(burgList[i], lineDir[i], nDir[i]);
      } else if (costheta2 <= (cosCritical2 * cosCritical2)) {
        /* mixed region */
        real8 acostheta2, temp;
        real8 nScrew[3], nNoScrew[3];
        VECTOR_COPY(nScrew, nDir[i]);
        NormalizedCrossVector(burgList[i], lineDir[i], nNoScrew);
        acostheta2 = acos(sqrt(costheta2));
        temp = cos(0.5 * (acostheta2 - angle1) * M_PI / (angle2 - angle1));
        temp = 1.0 - (temp * temp);
        for (j = 0; j < 3; j++) {
          nDir[i][j] = nScrew[j] + (nNoScrew[j] - nScrew[j]) * temp;
        }
      }
    }
  }

  /*
   *      It's possible this function was called for a node which had only zero-
   *      length segments (during SplitSurfaceNodes() for example).  If that is
   *      the case, just set the velocity to zero and return.
   */
  if (numNonZeroLenSegs == 0) {
    node->vX = 0.0;
    node->vY = 0.0;
    node->vZ = 0.0;
    return (0);
  }

  /*
   *      Now initialize the velocity
   */
  node->vX = 0.0e0;
  node->vY = 0.0e0;
  node->vZ = 0.0e0;

  /*
   *      Begin construction of the node drag matrix
   *
   *      Loop over all arms of the node, adding each arm's contribution
   *      to the drag matrix.
   */
  for (i = 0; i < numNbrs; i++) {
    real8 b[3], d[3], m[3], n[3], nCryst[3];

    /*
     *          If the segment is zero length (which can happen when
     *          the mobility function is being called from SplitMultiNodes())
     *          just skip the segment.
     */
    if (segLen[i] < eps) {
      continue;
    }

    VECTOR_COPY(b, burgList[i]);
    VECTOR_COPY(d, lineDir[i]);
    VECTOR_COPY(n, nDir[i]);

    halfLen = segLen[i] * 0.5;

    bDotb = DotProduct(b, b);
    costheta = DotProduct(d, b);
    costheta2 = costheta * costheta / bDotb;

    /*
     *          If needed, rotate a copy of the glide plane vector from the
     *          laboratory frame to the crystal frame.
     */
    if (param->useLabFrame) {
      Matrix33Vector3Multiply(home->rotMatrixInverse, n, nCryst);
    } else {
      VECTOR_COPY(nCryst, n);
    }

    /*
     *          arms not on [1 1 1] planes don't move as readily as
     *          other arms, so must be handled specially.
     */
    if (fabs(nCryst[X] * nCryst[Y] * nCryst[Z]) < 0.13608) {
      if (numNbrs == 2) {
        Btotal[0][0] += halfLen * Beclimb;
        Btotal[1][1] += halfLen * Beclimb;
        Btotal[2][2] += halfLen * Beclimb;
      } else {
        Btotal[0][0] += halfLen * (d[X] * d[X] * BlmBecl + Beclimb);
        Btotal[0][1] += halfLen * (d[X] * d[Y] * BlmBecl);
        Btotal[0][2] += halfLen * (d[X] * d[Z] * BlmBecl);
        Btotal[1][1] += halfLen * (d[Y] * d[Y] * BlmBecl + Beclimb);
        Btotal[1][2] += halfLen * (d[Y] * d[Z] * BlmBecl);
        Btotal[2][2] += halfLen * (d[Z] * d[Z] * BlmBecl + Beclimb);
      }
    } else {
      NormalizedCrossVector(n, d, m);
      Bglide = sqrt(invBedge2 + (invBscrew2 - invBedge2) * costheta2);
      Bglide = 1.0 / Bglide;
      Bclimb = Beclimb;
      Btotal[0][0] +=
          (m[X] * m[X] * Bglide + n[X] * n[X] * Bclimb + d[X] * d[X] * Bline) *
          halfLen;
      Btotal[0][1] +=
          (m[X] * m[Y] * Bglide + n[X] * n[Y] * Bclimb + d[X] * d[Y] * Bline) *
          halfLen;
      Btotal[0][2] +=
          (m[X] * m[Z] * Bglide + n[X] * n[Z] * Bclimb + d[X] * d[Z] * Bline) *
          halfLen;
      Btotal[1][1] +=
          (m[Y] * m[Y] * Bglide + n[Y] * n[Y] * Bclimb + d[Y] * d[Y] * Bline) *
          halfLen;
      Btotal[1][2] +=
          (m[Y] * m[Z] * Bglide + n[Y] * n[Z] * Bclimb + d[Y] * d[Z] * Bline) *
          halfLen;
      Btotal[2][2] +=
          (m[Z] * m[Z] * Bglide + n[Z] * n[Z] * Bclimb + d[Z] * d[Z] * Bline) *
          halfLen;
    }

  } /* End loop over arms */

  Btotal[1][0] = Btotal[0][1];
  Btotal[2][0] = Btotal[0][2];
  Btotal[2][1] = Btotal[1][2];

  /*
   *      At this point we should check if the matrix is invertable and
   *      if it isn't, find the eigen values and eigen vectors of the drag
   *      matrix, then invert the drag matrix keeping zero eigen values
   *      as zero.
   *
   *      FIX ME!  For now, we're assuming the matrix is invertible.
   */
  nForce[0] = node->fX;
  nForce[1] = node->fY;
  nForce[2] = node->fZ;

  if (Matrix33Invert(Btotal, invBtotal) == 0) {
    Fatal("Mobility_FCC_climb: Cannot invert 3X3 matrix!");
  }

  Matrix33Vector3Multiply(invBtotal, nForce, nVel);

  node->vX = nVel[0];
  node->vY = nVel[1];
  node->vZ = nVel[2];

  return (0);
}

int Mobility_FCC_0(Home_t *home, Node_t *node) {
  if (false) {
    return Mobility_FCC_climb(home, node);
  } else {
    return Mobility_FCC_Glide(home, node);
  }
}