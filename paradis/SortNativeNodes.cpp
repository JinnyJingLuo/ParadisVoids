/**************************************************************************
 *
 *  Function    : SortNativeNodes
 *  Description : Queue each native node onto the cell it currently falls into.
 *
 **************************************************************************/

#include "Home.h"
#include "Node.h"
#include "Cell.h"
#include "Util.h"

void AssignNodeToCell(Home_t *home, Node_t *node) {
  Param_t *param;
  int iCell, jCell, kCell, cellIdx;
  real8 probXmin, probYmin, probZmin;
  real8 cellXsize, cellYsize, cellZsize;
  Cell_t *cell;

  if (node == NULL) {
    return;
  }

  /*
   *      set the lower limit of the base area (excluding possible periodic
   * cells) and the size of each cell
   */
  param = home->param;
  probXmin = -0.5 * param->Dimensions[X];
  probYmin = -0.5 * param->Dimensions[Y];
  probZmin = -0.5 * param->Dimensions[Z];

  cellXsize = param->Dimensions[X] / param->nXcells;
  cellYsize = param->Dimensions[Y] / param->nYcells;
  cellZsize = param->Dimensions[Z] / param->nZcells;

  double dX = node->x;
  double dY = node->y;
  double dZ = node->z;

  GetInLocalCoordinates(home->param, dX, dY, dZ);
  /*
   *      put the node on its proper cell. If the index exceeds this domains
   *      range of native cells, put in the nearest native cell
   */
  iCell = (int)((dX - probXmin) / cellXsize);
  if (iCell < param->iCellNatMin)
    iCell = param->iCellNatMin;
  if (iCell > param->iCellNatMax)
    iCell = param->iCellNatMax;
  iCell++; /* compensate for periodic cells */

  jCell = (int)((dY - probYmin) / cellYsize);
  if (jCell < param->jCellNatMin)
    jCell = param->jCellNatMin;
  if (jCell > param->jCellNatMax)
    jCell = param->jCellNatMax;
  jCell++; /* compensate for periodic cells */

  kCell = (int)((dZ - probZmin) / cellZsize);
  if (kCell < param->kCellNatMin)
    kCell = param->kCellNatMin;
  if (kCell > param->kCellNatMax)
    kCell = param->kCellNatMax;
  kCell++; /* compensate for periodic cells */

  cellIdx = EncodeCellIdx(home, iCell, jCell, kCell);
  cell = home->cellKeys[cellIdx];
  node->nextInCell = cell->nodeQ;
  cell->nodeQ = node;
  cell->nodeCount++;

  node->cellIdx = cellIdx;
}

void SortNativeNodes(Home_t *home) {
  int i;
  Cell_t *cell;

  // initialize cell node lists
  for (i = 0; i < home->cellCount; i++) {
    cell = home->cellKeys[home->cellList[i]];
    cell->nodeQ = 0;
    cell->nodeCount = 0;
  }

  /* Loop thru active nodes, putting them in their proper cell. If the
   * index exceeds this domains range of native cells, put in the nearest
   * native cell
   */

  for (i = 0; i < home->newNodeKeyPtr; i++) {
    AssignNodeToCell(home, home->nodeKeys[i]);
  }
}
