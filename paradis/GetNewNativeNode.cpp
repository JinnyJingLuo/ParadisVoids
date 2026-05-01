/***************************************************************************
 *
 *      Function:    GetNewNativeNode
 *      Description: Pop a free Node_t from the free node Queue, and assign
 *                   it a free tag. NOTE: This routine does not allocate
 *                   any arrays specific to the number of neighbors. The
 *                   caller must use AllocNodeArms() to do this.
 *
 **************************************************************************/
#include "Home.h"
#include "Util.h"
#include "QueueOps.h"

Node_t *GetNewNativeNode(Home_t *home) {
  int newIdx;
  Node_t *newNode;

  newNode = PopFreeNodeQ(home);
  newIdx = GetFreeNodeTag(home);

  home->nodeKeys[newIdx] = newNode;

  newNode->myTag.domainID = home->myDomain;
  newNode->myTag.index = newIdx;
  newNode->cellIdx = -1;
  newNode->cell2Idx = -1;
  newNode->cell2QentIdx = -1;
  /*
   *      Explicitly zero out velocity so we don't end up
   *      with garbage when we are calculating the velocity
   *      delta between timesteps.
   */
  newNode->vX = 0.0;
  newNode->vY = 0.0;
  newNode->vZ = 0.0;

  newNode->dNx = 0.0;
  newNode->dNy = 0.0;
  newNode->dNz = 0.0;

  newNode->VelocityDapmingSteps = 0;

  newNode->flags = 0;
  // qjiao: clear interface flag;
  newNode->csState = ON_SLIP_PLANE;
  return (newNode);
}
