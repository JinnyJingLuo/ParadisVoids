/***************************************************************************
 *
 *      Module:       CellCharge
 *      Description:  Sum the contributions from each segment in a cell to
 *                    the total charge tensor in that cell, and distribute
 *                    the sum to all processors. The result is that all
 *                    processors have the net charge tensor of each cell in
 *                    the problem.
 *
 *      Includes functions:
 *
 *          CellCharge()
 *          FMSetTaylorExpansions()
 *          FMCellCharge()
 *          MonopoleCellCharge()
 *
 *      Last Modified: 04/08/2008 - Added explicit initialization of
 *                                  taylor coefficients for highest layer
 *                                  fmm cell if PBC is disabled.
 *
 ***************************************************************************/

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Home.h"
#include "Util.h"
#include "DSMPI.h"

void MonopoleCellCharge(Home_t *home) {

  int nCells, i, j, inode, inbr, iCell, jCell, kCell, cellIdx, chgIdx;
  real8 x1, y1, z1, dx, dy, dz;
  real8 b[3], dl[3];
  real8 *cellCharge;

  Param_t *param;
  Node_t *node, *nbr;
  Cell_t *cell;

  param = home->param;

  nCells = param->nXcells * param->nYcells * param->nZcells;

  /* Allocate and zero out a local charge array */

  cellCharge = (real8 *)malloc(9 * nCells * sizeof(real8));

  for (i = 0; i < 9 * nCells; i++) {
    cellCharge[i] = 0.0;
  }

  /* loop through the native nodes. For each node look at its neighbors.
   * If neighbor has a higher tag than node, then it is a non-redundant
   * segment.
   */

  for (inode = 0; inode < home->newNodeKeyPtr; inode++) {
    node = home->nodeKeys[inode];
    if (!node)
      continue;
    x1 = node->x;
    y1 = node->y;
    z1 = node->z;

    /* For consistency with the way segments are chosen for local force
     * calculation, a segment is considered to belong to the cell containing
     * the node owning the segment
     */

    cell = home->cellKeys[node->cellIdx];

    iCell = cell->xIndex;
    jCell = cell->yIndex;
    kCell = cell->zIndex;

    iCell--;
    jCell--;
    kCell--;

    cellIdx = kCell + param->nZcells * jCell +
              param->nZcells * param->nYcells * iCell;

    for (inbr = 0; inbr < node->numNbrs; inbr++) {

      nbr = GetNeighborNode(home, node, inbr);

      if (nbr == (Node_t *)NULL) {
        printf("WARNING: Neighbor not found at %s line %d\n", __FILE__,
               __LINE__);
        continue;
      }

      if (NodeOwnsSeg(home, node, nbr) == 0) {
        continue;
      }

      dx = nbr->x - x1;
      dy = nbr->y - y1;
      dz = nbr->z - z1;

      b[0] = node->burgX[inbr];
      b[1] = node->burgY[inbr];
      b[2] = node->burgZ[inbr];

      dl[0] = dx;
      dl[1] = dy;
      dl[2] = dz;

      /* Initialize chgIdx to the first element (charge[0][0]) of the tensor
       * for the segment's cell. The tensor is then looped thru in row order
       */

      chgIdx = 9 * cellIdx;
      for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
          cellCharge[chgIdx++] += b[i] * dl[j];
        }
      }

    } /* end for (inbr = 0 ; ...)  */
  }   /* end for (inode = 0 ; ...)  */

  /* Sum the individual cell charges over all processors (only a few at
   * most will contribute to any particular cell) and leave the sum for each
   * cell on all processors.
   */

  DSMPI::AllReduce(cellCharge, home->cellCharge, 9 * nCells, MPI_DOUBLE,
                   MPI_SUM);
  free(cellCharge);
}

void CellCharge(Home_t *home) {
  /*
   *      If we are not doing full n^2 force calculations, we need
   *      to prepare for the remote force calcs... but if we are
   *      doing full n^2 force calcs, this stuff is irrelevant.
   */
  MonopoleCellCharge(home);
  DSMPI::Barrier();
}
