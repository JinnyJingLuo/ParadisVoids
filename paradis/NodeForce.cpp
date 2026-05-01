/**************************************************************************
 *
 *      Module:  NodeForce -- This module contains various functions for
 *               calculating nodal forces
 *
 *      Includes functions:
 *               AddtoNodeForce()
 *               AddtoArmForce()
 *               ComputeFEMSegSigbRem()
 *               ComputeSegSigbRem()
 *               ExtPKForce()
 *               NodeForce()
 *               PKForce()
 *               SelfForce()
 *
 *
 *      NOTE:  This module contains several blocks of code which are
 *             conditionally compiled based on the _FEM and _FEMIMGSTRESS
 *             definitions.  These pre-processor definitions are only
 *             set when ParaDiS is being coupled with some locally developed
 *             finite element modules which are not included as part of
 *             the ParaDiS release.  Therefore, these blocks of code
 *             can be ignored by non-LLNL developers.
 *
 *      03/11/2009 - Replaced line tension model with new method from
 *                   Wei Cai, Sylvie Aubry, etc
 *
 *************************************************************************/
#include <stdio.h>
#include <math.h>

#ifdef PARALLEL
#include "mpi.h"
#endif

#include "Home.h"
#include "Util.h"
#include "Mobility.h"

#if defined _FEM | defined _FEMIMGSTRESS
#include "FEM.h"
#endif

/*
 *      Define a local structure in which to store some basic info we
 *      need when setting up a list of segment pairs for which forces
 *      need to be calculated
 */
typedef struct {
  Node_t *node1;
  Node_t *node2;
  Node_t *node3;
  Node_t *node4;
  real8 cellCenter[3];
} SegPair_t;

static real8 *cellCenterX = (real8 *)NULL;
static real8 *cellCenterY = (real8 *)NULL;
static real8 *cellCenterZ = (real8 *)NULL;

/*---------------------------------------------------------------------------
 *
 *    Function:        FreeCellCenters
 *    Description:     Release memory allocated to arrays for storing
 *                     the coordinates of the cell centers.  These arrays
 *                     are used in the GetFieldPointStressRem() function.
 *
 *-------------------------------------------------------------------------*/
void FreeCellCenters(void) {
  free(cellCenterX);
  free(cellCenterY);
  free(cellCenterZ);
}

/*---------------------------------------------------------------------------
 *
 *    Function:        AddtoNodeForce
 *    Description:     Increment the total nodal force for the indicated node
 *
 *-------------------------------------------------------------------------*/
void AddtoNodeForce(Node_t *node, real8 f[3]) {
  node->fX += f[0];
  node->fY += f[1];
  node->fZ += f[2];
}

/*---------------------------------------------------------------------------
 *
 *    Function:        AddtoArmForce
 *    Description:     Increment the nodal force attributable to
 *                     a specific arm of a node by the specified
 *                     amount.
 *
 *-------------------------------------------------------------------------*/
void AddtoArmForce(Node_t *node, int arm, real8 f[3]) {
  node->armfx[arm] += f[0];
  node->armfy[arm] += f[1];
  node->armfz[arm] += f[2];
}

void ZeroNodeForces(Home_t *home) {
  int i, j;
  Node_t *node, *nbr;

  for (i = 0; i < home->newNodeKeyPtr; i++) {
    if ((node = home->nodeKeys[i]) == NULL) {
      continue;
    }

    node->fX = 0.0;
    node->fY = 0.0;
    node->fZ = 0.0;

    for (j = 0; j < node->numNbrs; j++) {
      if ((nbr = GetNeighborNode(home, node, j)) == NULL) {
        continue;
      }
      node->armfx[j] = 0.0;
      node->armfy[j] = 0.0;
      node->armfz[j] = 0.0;
      node->sigbRem[3 * j] = 0.0;
      node->sigbRem[3 * j + 1] = 0.0;
      node->sigbRem[3 * j + 2] = 0.0;
    }
  }

  /*
   *      Now do the same for the ghost nodes.
   */
  node = home->ghostNodeQ;
  while (node) {
    node->fX = 0.0;
    node->fY = 0.0;
    node->fZ = 0.0;

    for (j = 0; j < node->numNbrs; j++) {
      if ((nbr = GetNeighborNode(home, node, j)) == NULL) {
        continue;
      }
      node->armfx[j] = 0.0;
      node->armfy[j] = 0.0;
      node->armfz[j] = 0.0;
      node->sigbRem[3 * j] = 0.0;
      node->sigbRem[3 * j + 1] = 0.0;
      node->sigbRem[3 * j + 2] = 0.0;
    }
    node = node->next;
  }
}

static void GetFieldPointStressRem(Home_t *home, real8 x, real8 y, real8 z,
                                   int cellX, int cellY, int cellZ,
                                   real8 totStress[3][3]) {
  int i, j, includePrimary;
  int cx, cy, cz, cellIndex;
  int xSkip1, xSkip2, xSkip3;
  int ySkip1, ySkip2, ySkip3;
  int zSkip1, zSkip2, zSkip3;
  real8 dx, dy, dz;
  real8 xc, yc, zc;
  real8 burgX, burgY, burgZ;
  real8 delSig[3][3];
  Param_t *param;

  param = home->param;

  /*
   *      First time into this function, allocate and initialize some
   *      static arrays
   */
  if (cellCenterX == NULL) {
    real8 Lx, Ly, Lz;
    real8 cellXsize, cellYsize, cellZsize, xStart, yStart, zStart;

    Lx = param->Lx;
    Ly = param->Ly;
    Lz = param->Lz;

    cellCenterX = (real8 *)malloc(param->nXcells * sizeof(real8));
    cellCenterY = (real8 *)malloc(param->nYcells * sizeof(real8));
    cellCenterZ = (real8 *)malloc(param->nZcells * sizeof(real8));

    cellXsize = Lx / param->nXcells;
    cellYsize = Ly / param->nYcells;
    cellZsize = Lz / param->nZcells;

    xStart = 0.5 * (-param->Dimensions[X] + cellXsize);
    yStart = 0.5 * (-param->Dimensions[Y] + cellYsize);
    zStart = 0.5 * (-param->Dimensions[Z] + cellZsize);

    for (i = 0; i < param->nXcells; i++) {
      cellCenterX[i] = xStart + i * cellXsize;
    }
    for (i = 0; i < param->nYcells; i++) {
      cellCenterY[i] = yStart + i * cellYsize;
    }
    for (i = 0; i < param->nZcells; i++) {
      cellCenterZ[i] = zStart + i * cellZsize;
    }
  }

  /*
   *      Initialize the stress at the field point
   */
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      totStress[i][j] = 0.0;
    }
  }

  /*
   *      Since we have to treat neighboring and remote cells
   *      differently, we need to identify which cells are which...
   *      If the neighbor cells are outside the primary image, wrap
   *      the cell indices around to the other side of the box ONLY
   *      if periodic boundaries are enabld!
   */
  xSkip1 = cellX - 1;
  if (xSkip1 < 0) {
    xSkip1 = 0;
  }
  xSkip2 = cellX;
  xSkip3 = cellX + 1;
  if (xSkip3 >= param->nXcells) {
    xSkip3 = param->nXcells - 1;
  }

  ySkip1 = cellY - 1;
  if (ySkip1 < 0) {
    ySkip1 = 0;
  }
  ySkip2 = cellY;
  ySkip3 = cellY + 1;
  if (ySkip3 >= param->nYcells) {
    ySkip3 = param->nYcells - 1;
  }

  zSkip1 = cellZ - 1;
  if (zSkip1 < 0) {
    zSkip1 = 0;
  }
  zSkip2 = cellZ;
  zSkip3 = cellZ + 1;
  if (zSkip3 >= param->nZcells) {
    zSkip3 = param->nZcells - 1;
  }

  /*
   *      Loop through all cells, and add the stress contribution from
   *      the cell to stress at the segment midpoint.
   */
  for (cx = 0; cx < param->nXcells; cx++) {
    for (cy = 0; cy < param->nYcells; cy++) {
      for (cz = 0; cz < param->nZcells; cz++) {
        includePrimary = !((cx == xSkip1 || cx == xSkip2 || cx == xSkip3) &&
                           (cy == ySkip1 || cy == ySkip2 || cy == ySkip3) &&
                           (cz == zSkip1 || cz == zSkip2 || cz == zSkip3));

        /*
         *              Get the center point of cell [cx, cy, cz]
         */
        xc = cellCenterX[cx];
        yc = cellCenterY[cy];
        zc = cellCenterZ[cz];

        /*
         *              Get the stress at the specified point caused
         *              by the net charge tensor of the current cell.
         */
        dx = xc - x;
        dy = yc - y;
        dz = zc - z;

        GetPrimaryImage(param, &dx, &dy, &dz);

        xc = x + dx;
        yc = y + dy;
        zc = z + dz;

        cellIndex =
            cz + param->nZcells * cy + param->nZcells * param->nYcells * cx;

        /*
         *              Stress (charge[.,1], [1,0,0])
         */
        burgX = home->cellCharge[9 * cellIndex];
        burgY = home->cellCharge[9 * cellIndex + 3];
        burgZ = home->cellCharge[9 * cellIndex + 6];

        dx = 1;
        dy = 0;
        dz = 0;

        dSegImgStress(home, delSig, xc, yc, zc, dx, dy, dz, burgX, burgY, burgZ,
                      x, y, z, includePrimary);

        for (i = 0; i < 3; i++) {
          for (j = 0; j < 3; j++) {
            totStress[i][j] += delSig[i][j];
          }
        }
        /*
         *              Stress (charge[.,2], [0,1,0])
         */
        burgX = home->cellCharge[9 * cellIndex + 1];
        burgY = home->cellCharge[9 * cellIndex + 4];
        burgZ = home->cellCharge[9 * cellIndex + 7];

        dx = 0;
        dy = 1;
        dz = 0;

        dSegImgStress(home, delSig, xc, yc, zc, dx, dy, dz, burgX, burgY, burgZ,
                      x, y, z, includePrimary);

        for (i = 0; i < 3; i++) {
          for (j = 0; j < 3; j++) {
            totStress[i][j] += delSig[i][j];
          }
        }

        /*
         *              Stress (charge[.,3], [0,0,1])
         */
        burgX = home->cellCharge[9 * cellIndex + 2];
        burgY = home->cellCharge[9 * cellIndex + 5];
        burgZ = home->cellCharge[9 * cellIndex + 8];

        dx = 0;
        dy = 0;
        dz = 1;

        dSegImgStress(home, delSig, xc, yc, zc, dx, dy, dz, burgX, burgY, burgZ,
                      x, y, z, includePrimary);

        for (i = 0; i < 3; i++) {
          for (j = 0; j < 3; j++) {
            totStress[i][j] += delSig[i][j];
          }
        }
      } /* end for(cz = 0; ...) */
    }   /* end for(cy = 0; ...) */
  }     /* end for(cx = 0; ...) */
}

/*-------------------------------------------------------------------------
 *
 *      Function:     ComputeSigSegbRem
 *      Description:  For each native segment, calculate the resultant
 *                    stress at the segment midpoint due to the cellCharges
 *                    (accumulated bXdl in a cell) from all cells. Store
 *                    the dot product of this tensor with the segment's
 *                    burger's vector (i.e sigma dot b) in the segment's
 *                    remote sigb.
 *
 *                    NOTE: For remote cells, the resultant stress on each
 *                          segment includes contributions from both the
 *                          primary and periodic images.  For the segment's
 *                          local and neighboring cells, this stress includes
 *                          contributions only from the periodic images.
 *                          Force from interactions with the primary image
 *                          have been computed elsewhere.
 *      Arguments:
 *          reqType   Indicates the type of calculation being requested.
 *                    Permitted values are:
 *
 *                        PARTIAL  Remote sigb is calculated only for nodes
 *                                 flagged to have forces recalculated.
 *                        FULL     Remote sigb is calculated for all nodes.
 *
 *-----------------------------------------------------------------------*/
void ComputeSegSigbRem(Home_t *home) {
  int inode, ti, myXcell, myYcell, myZcell;
  int armID1, armID2;
  real8 x1, y1, z1, dx, dy, dz, xm, ym, zm;
  real8 bx1, by1, bz1, sigb1, sigb2, sigb3;
  real8 totRemSig[3][3];
  Param_t *param;
  Node_t *node, *nbr;
  Cell_t *cell;

  param = home->param;

  /*
   *      loop thru the native nodes
   */
  for (inode = 0; inode < home->newNodeKeyPtr; inode++) {
    if ((node = home->nodeKeys[inode]) == NULL) {
      continue;
    }
    /*
     *          loop thru the native node's arms.  If the segment is owned
     *          by the neighboring node, don't do the calculation in this
     *          iteration of the loop.
     */
    for (ti = 0; ti < node->numNbrs; ti++) {
      nbr = GetNeighborNode(home, node, ti);
      if (nbr == NULL) {
        continue;
      }

      if (NodeOwnsSeg(home, node, nbr) == 0) {
        continue;
      }

      armID1 = GetArmID(node, nbr);
      armID2 = GetArmID(nbr, node);

      x1 = node->x;
      y1 = node->y;
      z1 = node->z;

      /*
       *              Initialize the sigb on the arm
       */
      node->sigbRem[3 * armID1] = 0.0;
      node->sigbRem[3 * armID1 + 1] = 0.0;
      node->sigbRem[3 * armID1 + 2] = 0.0;

      /*
       *              Get the midpoint of the segment
       */
      dx = nbr->x - x1;
      dy = nbr->y - y1;
      dz = nbr->z - z1;

      GetPrimaryImage(param, &dx, &dy, &dz);

      xm = x1 + dx * 0.5;
      ym = y1 + dy * 0.5;
      zm = z1 + dz * 0.5;

      /*
       *              Get the cell indices for the cell containing the node
       *              and remove adjustments for periodic cells, then find
       *              the stress at the specied point from all segments in
       *              remote cells.
       */
      cell = home->cellKeys[node->cellIdx];

      myXcell = cell->xIndex;
      myYcell = cell->yIndex;
      myZcell = cell->zIndex;

      myXcell--;
      myYcell--;
      myZcell--;

      GetFieldPointStressRem(home, xm, ym, zm, myXcell, myYcell, myZcell,
                             totRemSig);

      /*
       *              Calculate Sigma dot b from remote stess on this segment,
       * and store in the node arm
       */
      bx1 = node->burgX[armID1];
      by1 = node->burgY[armID1];
      bz1 = node->burgZ[armID1];

      sigb1 =
          totRemSig[0][0] * bx1 + totRemSig[0][1] * by1 + totRemSig[0][2] * bz1;
      sigb2 =
          totRemSig[1][0] * bx1 + totRemSig[1][1] * by1 + totRemSig[1][2] * bz1;
      sigb3 =
          totRemSig[2][0] * bx1 + totRemSig[2][1] * by1 + totRemSig[2][2] * bz1;

      node->sigbRem[3 * armID1] = sigb1;
      node->sigbRem[3 * armID1 + 1] = sigb2;
      node->sigbRem[3 * armID1 + 2] = sigb3;

      nbr->sigbRem[3 * armID2] = sigb1;
      nbr->sigbRem[3 * armID2 + 1] = sigb2;
      nbr->sigbRem[3 * armID2 + 2] = sigb3;
    } /* for(ti=0;ti<nc;ti++) */
  }   /* for(inode=0;...) */
}

/*-------------------------------------------------------------------------
 *
 *      Function:     NodeForce
 *      Description:  This function does no force calulcations directly
 *                    but drives the calculations via lower level functions.
 *                    It generates the sig.b from remote segments, the
 *                    forces from direct segment-to-segment interactions,
 *                    then for each local segment, calculates the self-force,
 *                    the force from extrenal stress, and the force from
 *                    the remote sig.b.
 *
 *-----------------------------------------------------------------------*/
void NodeForce(Home_t *home, ParadisExternalLoadServer *poExternalLoadServer,
               ParadisPrecipitateServer *poPrecipitateServer) {
  int i, nc, ti, nbrArm;
  real8 f1[3], f2[3];
  Node_t *node, *nbr;
  Param_t *param;

  param = home->param;
  /*
   *      Reset all node forces to zero (local and ghost nodes)
   */
  ZeroNodeForces(home);
/*
 *      NOTE: If full n^2 seg/seg forces are being calculated, we don't
 *      do any remote force calcs
 */

#ifndef FULL_N2_FORCES
  ComputeSegSigbRem(home);
#endif

  /*
   *      Now handle all the force calculations that must be done by the
   *      local domain.  This includes self-force, PK force, and far-field
   *      forces for all native segments plus any segment-to-segment
   *      interactions for segment pairs 'owned' by this domain.
   *
   *      All force calulations for a given segment will be calculated
   *      only once, and the calculating domain will distribute calculated
   *      forces to remote domains as necessary.
   */
  LocalSegForces(home, poExternalLoadServer, poPrecipitateServer);

#if PARALLEL
  MPI_Barrier(MPI_COMM_WORLD);
#endif
}

void PKForce(real8 sigb[3], real8 x1, real8 y1, real8 z1, real8 x2, real8 y2,
             real8 z2, real8 f1[3], real8 f2[3]) {
  real8 xi[3], ft[3];
  int j;

  xi[0] = x2 - x1;
  xi[1] = y2 - y1;
  xi[2] = z2 - z1;

  cross(sigb, xi, ft);

  for (j = 0; j < 3; j++) {
    f1[j] = ft[j] * 0.5;
    f2[j] = ft[j] * 0.5;
  }
}

/*-------------------------------------------------------------------------
 *
 *      Function:       SelfForce
 *      Description:    only the isotropic expression for self forces is
 *supportedhere
 *
 *      Arguments:
 *
 *------------------------------------------------------------------------*/
void SelfForce(int coreOnly, real8 MU, real8 NU, real8 bx, real8 by, real8 bz,
               real8 x1, real8 y1, real8 z1, real8 x2, real8 y2, real8 z2,
               real8 a, real8 Ecore, real8 f1[3], real8 f2[3]) {

  f1[0] = 0.0;
  f1[1] = 0.0;
  f1[2] = 0.0;

  f2[0] = 0.0;
  f2[1] = 0.0;
  f2[2] = 0.0;

  real8 tx, ty, tz, L, La, S;
  real8 bs, bs2, bex, bey, bez, be2, fL, ft;

  GetUnitVector(1, x1, y1, z1, x2, y2, z2, &tx, &ty, &tz, &L);
  if (L < 1.0E-6) {
    return;
  }
  bs = bx * tx + by * ty + bz * tz;
  bex = bx - bs * tx;
  bey = by - bs * ty;
  bez = bz - bs * tz;
  be2 = (bex * bex + bey * bey + bez * bez);
  bs2 = bs * bs;

  La = sqrt(L * L + a * a);

  if (coreOnly) {
    S = 0.0;
  } else {
    S = (-(2 * NU * La + (1 - NU) * a * a / La - (1 + NU) * a) / L +
         (NU * log((La + L) / a) - (1 - NU) * 0.5 * L / La)) *
        MU / 4 / M_PI / (1 - NU) * bs;
  }

  /* Ecore = MU/(4*pi) log(a/a0) */
  /* Account for the self force due to the core energy change when
     the core radius changes from a0-->a, M. Tang, 7/19/2004 */
  fL = -Ecore * (bs2 + be2 / (1 - NU));
  ft = Ecore * 2 * bs * NU / (1 - NU);

  f2[0] = bex * (S + ft) + fL * tx;
  f2[1] = bey * (S + ft) + fL * ty;
  f2[2] = bez * (S + ft) + fL * tz;

  f1[0] = -f2[0];
  f1[1] = -f2[1];
  f1[2] = -f2[2];
}

/**************************************************************************
 *
 *      Function:    StressDueToSeg
 *      Description: Calculate the stress at point p from the segment
 *                   starting at point p1 and ending at point p2.
 *
 *      Arguments:
 *         px, py, pz     coordinates of field point at which stress is to
 *                        be evaluated
 *         p1x, p1y, p1z  starting position of the dislocation segment
 *         p2x, p2y, p2z  ending position of the dislocation segment
 *         bx, by, bz     burgers vector associated with segment going
 *                        from p1 to p2
 *         a              core value
 *         MU             shear modulus
 *         NU             poisson ratio
 *         stress         array of stresses form the indicated segment
 *                        at the field point requested
 *                            [0] = stressxx
 *                            [1] = stressyy
 *                            [2] = stresszz
 *                            [3] = stressxy
 *                            [4] = stressyz
 *                            [5] = stressxz
 *
 *************************************************************************/
void StressDueToSeg(real8 px, real8 py, real8 pz, real8 p1x, real8 p1y,
                    real8 p1z, real8 p2x, real8 p2y, real8 p2z, real8 bx,
                    real8 by, real8 bz, real8 a, real8 MU, real8 NU,
                    real8 *stress) {
  real8 oneoverLp, common;
  real8 vec1x, vec1y, vec1z;
  real8 tpx, tpy, tpz;
  real8 Rx, Ry, Rz, Rdt;
  real8 ndx, ndy, ndz;
  real8 d2, s1, s2, a2, a2_d2, a2d2inv;
  real8 Ra, Rainv, Ra3inv, sRa3inv;
  real8 s_03a, s_13a, s_05a, s_15a, s_25a;
  real8 s_03b, s_13b, s_05b, s_15b, s_25b;
  real8 s_03, s_13, s_05, s_15, s_25;
  real8 m4p, m8p, m4pn, mn4pn, a2m8p;
  real8 txbx, txby, txbz;
  real8 dxbx, dxby, dxbz;
  real8 dxbdt, dmdxx, dmdyy, dmdzz, dmdxy, dmdyz, dmdxz;
  real8 tmtxx, tmtyy, tmtzz, tmtxy, tmtyz, tmtxz;
  real8 tmdxx, tmdyy, tmdzz, tmdxy, tmdyz, tmdxz;
  real8 tmtxbxx, tmtxbyy, tmtxbzz, tmtxbxy, tmtxbyz, tmtxbxz;
  real8 dmtxbxx, dmtxbyy, dmtxbzz, dmtxbxy, dmtxbyz, dmtxbxz;
  real8 tmdxbxx, tmdxbyy, tmdxbzz, tmdxbxy, tmdxbyz, tmdxbxz;
  real8 I_03xx, I_03yy, I_03zz, I_03xy, I_03yz, I_03xz;
  real8 I_13xx, I_13yy, I_13zz, I_13xy, I_13yz, I_13xz;
  real8 I_05xx, I_05yy, I_05zz, I_05xy, I_05yz, I_05xz;
  real8 I_15xx, I_15yy, I_15zz, I_15xy, I_15yz, I_15xz;
  real8 I_25xx, I_25yy, I_25zz, I_25xy, I_25yz, I_25xz;

  vec1x = p2x - p1x;
  vec1y = p2y - p1y;
  vec1z = p2z - p1z;

  oneoverLp = 1 / sqrt(vec1x * vec1x + vec1y * vec1y + vec1z * vec1z);

  tpx = vec1x * oneoverLp;
  tpy = vec1y * oneoverLp;
  tpz = vec1z * oneoverLp;

  Rx = px - p1x;
  Ry = py - p1y;
  Rz = pz - p1z;

  Rdt = Rx * tpx + Ry * tpy + Rz * tpz;

  ndx = Rx - Rdt * tpx;
  ndy = Ry - Rdt * tpy;
  ndz = Rz - Rdt * tpz;

  d2 = ndx * ndx + ndy * ndy + ndz * ndz;

  s1 = -Rdt;
  s2 = -((px - p2x) * tpx + (py - p2y) * tpy + (pz - p2z) * tpz);
  a2 = a * a;
  a2_d2 = a2 + d2;
  a2d2inv = 1 / a2_d2;

  Ra = sqrt(a2_d2 + s1 * s1);
  Rainv = 1 / Ra;
  Ra3inv = Rainv * Rainv * Rainv;
  sRa3inv = s1 * Ra3inv;

  s_03a = s1 * Rainv * a2d2inv;
  s_13a = -Rainv;
  s_05a = (2 * s_03a + sRa3inv) * a2d2inv;
  s_15a = -Ra3inv;
  s_25a = s_03a - sRa3inv;

  Ra = sqrt(a2_d2 + s2 * s2);
  Rainv = 1 / Ra;
  Ra3inv = Rainv * Rainv * Rainv;
  sRa3inv = s2 * Ra3inv;

  s_03b = s2 * Rainv * a2d2inv;
  s_13b = -Rainv;
  s_05b = (2 * s_03b + sRa3inv) * a2d2inv;
  s_15b = -Ra3inv;
  s_25b = s_03b - sRa3inv;

  s_03 = s_03b - s_03a;
  s_13 = s_13b - s_13a;
  s_05 = s_05b - s_05a;
  s_15 = s_15b - s_15a;
  s_25 = s_25b - s_25a;

  m4p = 0.25 * MU / M_PI;
  m8p = 0.5 * m4p;
  m4pn = m4p / (1 - NU);
  mn4pn = m4pn * NU;
  a2m8p = a2 * m8p;

  txbx = tpy * bz - tpz * by;
  txby = tpz * bx - tpx * bz;
  txbz = tpx * by - tpy * bx;

  dxbx = ndy * bz - ndz * by;
  dxby = ndz * bx - ndx * bz;
  dxbz = ndx * by - ndy * bx;

  dxbdt = dxbx * tpx + dxby * tpy + dxbz * tpz;

  dmdxx = ndx * ndx;
  dmdyy = ndy * ndy;
  dmdzz = ndz * ndz;
  dmdxy = ndx * ndy;
  dmdyz = ndy * ndz;
  dmdxz = ndx * ndz;

  tmtxx = tpx * tpx;
  tmtyy = tpy * tpy;
  tmtzz = tpz * tpz;
  tmtxy = tpx * tpy;
  tmtyz = tpy * tpz;
  tmtxz = tpx * tpz;

  tmdxx = 2 * tpx * ndx;
  tmdyy = 2 * tpy * ndy;
  tmdzz = 2 * tpz * ndz;
  tmdxy = tpx * ndy + tpy * ndx;
  tmdyz = tpy * ndz + tpz * ndy;
  tmdxz = tpx * ndz + tpz * ndx;

  tmtxbxx = 2 * tpx * txbx;
  tmtxbyy = 2 * tpy * txby;
  tmtxbzz = 2 * tpz * txbz;
  tmtxbxy = tpx * txby + tpy * txbx;
  tmtxbyz = tpy * txbz + tpz * txby;
  tmtxbxz = tpx * txbz + tpz * txbx;

  dmtxbxx = 2 * ndx * txbx;
  dmtxbyy = 2 * ndy * txby;
  dmtxbzz = 2 * ndz * txbz;
  dmtxbxy = ndx * txby + ndy * txbx;
  dmtxbyz = ndy * txbz + ndz * txby;
  dmtxbxz = ndx * txbz + ndz * txbx;

  tmdxbxx = 2 * tpx * dxbx;
  tmdxbyy = 2 * tpy * dxby;
  tmdxbzz = 2 * tpz * dxbz;
  tmdxbxy = tpx * dxby + tpy * dxbx;
  tmdxbyz = tpy * dxbz + tpz * dxby;
  tmdxbxz = tpx * dxbz + tpz * dxbx;

  common = m4pn * dxbdt;

  I_03xx = common + m4pn * dmtxbxx - m4p * tmdxbxx;
  I_03yy = common + m4pn * dmtxbyy - m4p * tmdxbyy;
  I_03zz = common + m4pn * dmtxbzz - m4p * tmdxbzz;
  I_03xy = m4pn * dmtxbxy - m4p * tmdxbxy;
  I_03yz = m4pn * dmtxbyz - m4p * tmdxbyz;
  I_03xz = m4pn * dmtxbxz - m4p * tmdxbxz;

  I_13xx = -mn4pn * tmtxbxx;
  I_13yy = -mn4pn * tmtxbyy;
  I_13zz = -mn4pn * tmtxbzz;
  I_13xy = -mn4pn * tmtxbxy;
  I_13yz = -mn4pn * tmtxbyz;
  I_13xz = -mn4pn * tmtxbxz;

  I_05xx = common * (a2 + dmdxx) - a2m8p * tmdxbxx;
  I_05yy = common * (a2 + dmdyy) - a2m8p * tmdxbyy;
  I_05zz = common * (a2 + dmdzz) - a2m8p * tmdxbzz;
  I_05xy = common * dmdxy - a2m8p * tmdxbxy;
  I_05yz = common * dmdyz - a2m8p * tmdxbyz;
  I_05xz = common * dmdxz - a2m8p * tmdxbxz;

  I_15xx = a2m8p * tmtxbxx - common * tmdxx;
  I_15yy = a2m8p * tmtxbyy - common * tmdyy;
  I_15zz = a2m8p * tmtxbzz - common * tmdzz;
  I_15xy = a2m8p * tmtxbxy - common * tmdxy;
  I_15yz = a2m8p * tmtxbyz - common * tmdyz;
  I_15xz = a2m8p * tmtxbxz - common * tmdxz;

  I_25xx = common * tmtxx;
  I_25yy = common * tmtyy;
  I_25zz = common * tmtzz;
  I_25xy = common * tmtxy;
  I_25yz = common * tmtyz;
  I_25xz = common * tmtxz;

  stress[0] = I_03xx * s_03 + I_13xx * s_13 + I_05xx * s_05 + I_15xx * s_15 +
              I_25xx * s_25;

  stress[1] = I_03yy * s_03 + I_13yy * s_13 + I_05yy * s_05 + I_15yy * s_15 +
              I_25yy * s_25;

  stress[2] = I_03zz * s_03 + I_13zz * s_13 + I_05zz * s_05 + I_15zz * s_15 +
              I_25zz * s_25;

  stress[3] = I_03xy * s_03 + I_13xy * s_13 + I_05xy * s_05 + I_15xy * s_15 +
              I_25xy * s_25;

  stress[4] = I_03yz * s_03 + I_13yz * s_13 + I_05yz * s_05 + I_15yz * s_15 +
              I_25yz * s_25;

  stress[5] = I_03xz * s_03 + I_13xz * s_13 + I_05xz * s_05 + I_15xz * s_15 +
              I_25xz * s_25;

  return;
}

/*-------------------------------------------------------------------------
 *
 *      Function:       GetFieldPointStress
 *      Description:
 *
 *      Arguments:
 *
 *------------------------------------------------------------------------*/
void GetFieldPointStress(Home_t *home, real8 x, real8 y, real8 z,
                         real8 totStress[3][3]) {
  int cellX, cellY, cellZ, cellIndex, arm;
  int cx, cy, cz;
  int minXIndex, minYIndex, minZIndex;
  int maxXIndex, maxYIndex, maxZIndex;
  real8 Lx, Ly, Lz;
  real8 xc, yc, zc;
  real8 a, MU, NU;
  real8 bx, by, bz;
  real8 p1x, p1y, p1z;
  real8 p2x, p2y, p2z;
  real8 stress[6];
  Node_t *node1, *node2;
  Cell_t *cell;
  Param_t *param;

  param = home->param;

  a = param->rc;
  MU = param->shearModulus;
  NU = param->pois;

  totStress[0][0] = 0.0;
  totStress[0][1] = 0.0;
  totStress[0][2] = 0.0;
  totStress[1][0] = 0.0;
  totStress[1][1] = 0.0;
  totStress[1][2] = 0.0;
  totStress[2][0] = 0.0;
  totStress[2][1] = 0.0;
  totStress[2][2] = 0.0;

  /*
   *      Get the indices of the cell containing the field point and
   *      determine the block of cells immediately neighboring that
   *      cell.
   */
  Lx = param->Lx;
  Ly = param->Ly;
  Lz = param->Lz;

  cellX = (x + 0.5 * param->Dimensions[X]) / (Lx / param->nXcells) + 1;
  cellY = (y + 0.5 * param->Dimensions[Y]) / (Ly / param->nYcells) + 1;
  cellZ = (z + 0.5 * param->Dimensions[Z]) / (Lz / param->nZcells) + 1;

  /*
   *      Determine the minimum and maximum cell indices (in each
   *      dimension) of the block of cells encompassing the cell
   *      containing the field point, and all immediate neighbors
   *      of that cell.  Allow for periodic boundary conditions.
   */

  minXIndex = MAX(1, cellX - 1);
  maxXIndex = MIN(param->nXcells, cellX + 1);

  minYIndex = MAX(1, cellY - 1);
  maxYIndex = MIN(param->nYcells, cellY + 1);

  minZIndex = MAX(1, cellZ - 1);
  maxZIndex = MIN(param->nZcells, cellZ + 1);

  /*
   *      Loop though all the cells in the block.
   */
  for (cx = minXIndex; cx <= maxXIndex; cx++) {
    for (cy = minYIndex; cy <= maxYIndex; cy++) {
      for (cz = minZIndex; cz <= maxZIndex; cz++) {

        cellIndex = EncodeCellIdx(home, cx, cy, cz);
        cell = home->cellKeys[cellIndex];

        if (cell == (Cell_t *)NULL)
          continue;

        /*
         *                  Find the center of this cell and convert it to the
         *                  point in the image closest to the field point.
         *                  We need use the cell center as the base because
         *                  all the Burgers vectors of segments inside one cell
         *                  (belonging to the same PBC image) are summed
         *                  together for 'remote' stress calculations, so
         *                  cannot separate them for 'local' stress
         * calculations.
         */
        FindCellCenter(param, (real8)(cx - 1), (real8)(cy - 1), (real8)(cz - 1),
                       2, &xc, &yc, &zc);

        /*
         *                  Loop over all nodes in this cell and over each
         * segment attached to the node.  Skip any segment that is not owned by
         * node1.
         */
        node1 = cell->nodeQ;

        for (; node1 != (Node_t *)NULL; node1 = node1->nextInCell) {
          for (arm = 0; arm < node1->numNbrs; arm++) {

            node2 = GetNeighborNode(home, node1, arm);

            if (node2 == (Node_t *)NULL)
              continue;
            if (OrderNodes(node1, node2) >= 0)
              continue;

            p1x = node1->x;
            p1y = node1->y;
            p1z = node1->z;

            p2x = node2->x;
            p2y = node2->y;
            p2z = node2->z;

            bx = node1->burgX[arm];
            by = node1->burgY[arm];
            bz = node1->burgZ[arm];

            StressDueToSeg(x, y, z, p1x, p1y, p1z, p2x, p2y, p2z, bx, by, bz, a,
                           MU, NU, stress);

            totStress[0][0] += stress[0];
            totStress[1][1] += stress[1];
            totStress[2][2] += stress[2];
            totStress[0][1] += stress[3];
            totStress[1][2] += stress[4];
            totStress[0][2] += stress[5];
          }
        }
      }
    }
  }

  totStress[1][0] = totStress[0][1];
  totStress[2][0] = totStress[0][2];
  totStress[2][1] = totStress[1][2];
}
