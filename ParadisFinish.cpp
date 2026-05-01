/*-------------------------------------------------------------------------
 *
 *      Function:     ParadisFinish
 *      Description:  Handles pre-termination processing including
 *                    output generation, cleanup of any X-Window
 *                    display, release of any remaining dynamically
 *                    allocated memory, etc.
 *
 *-----------------------------------------------------------------------*/
#include "Home.h"
#include "Init.h"
#include "Util.h"
#include "Decomp.h"
#include "DSMPI.h"

/**************************************************************************
 *
 *      Function:     ReleaseMemory
 *      Description:  This function is called during the program
 *                    termination procedure to release all remaining
 *                    dynamically allocated memory.  This is not truly
 *                    necessary, however, doing so facilitates the
 *                    process of hunting memory leaks by cleaning up
 *                    all the known allocated memory blocks.
 *
 *************************************************************************/
void ReleaseMemory(Home_t *home) {
  int i, j, k, domID;
  Node_t *node;
  Param_t *param;
  NodeBlock_t *nodeBlock, *thisNodeBlock;

  if (home == (Home_t *)NULL) {
    return;
  }

  param = home->param;

  /*
   *      Use existing functions (if available) to clean up
   */

  DLBfreeOld(home);
  FreeRijm();
  FreeRijmPBC();

  /*
   *      Loop through all allocated blocks of node structures.  Free
   *      any arrays associated with individual nodes then free up the
   *      blocks of node structures.
   */
  nodeBlock = home->nodeBlockQ;

  while (nodeBlock) {

    node = nodeBlock->nodes;

    for (i = 0; i < NODE_BLOCK_COUNT; i++) {

      if (node->numNbrs) {
        FreeNodeArms(node);
      }
      node++;
    }

    thisNodeBlock = nodeBlock;
    nodeBlock = nodeBlock->next;

    memset(thisNodeBlock->nodes, 0, NODE_BLOCK_COUNT * sizeof(Node_t));
    free(thisNodeBlock->nodes);
    free(thisNodeBlock);
  }

  home->nodeBlockQ = 0;

  if (home->nodeKeys) {
    free(home->nodeKeys);
    home->nodeKeys = NULL;
  }
  home->newNodeKeyPtr = 0;
  home->newNodeKeyMax = 0;
  home->nativeNodeQ = NULL;
  home->ghostNodeQ = NULL;
  home->freeNodeQ = NULL;
  home->lastFreeNode = NULL;
  home->lastGhostNode = NULL;
  home->nodeBlockQ = NULL;

  if (home->param->SType > 0) {
    if (home->stress) {
      for (i = 0; i < param->GP_x; i++) {
        for (j = 0; j < param->GP_y; j++) {
          for (k = 0; k < param->GP_z; k++) {
            free(home->stress[i][j][k]);
          }
          free(home->stress[i][j]);
        }
        free(home->stress[i]);
      }
      free(home->stress);
      home->stress = NULL;
    }
  }

  /*
   *      Free the heap used for recycling nodes that have bee deleted
   */
  if (home->recycledNodeHeapSize > 0) {
    free(home->recycledNodeHeap);
    home->recycledNodeHeap = (int *)NULL;
  }
  home->recycledNodeHeapSize = 0;
  home->recycledNodeHeapEnts = 0;

  /*
   *      Free any memory associated with the list of topological
   *      operations distributed to remote domains.
   */
  if (home->opList) {
    free(home->opList);
    home->opList = NULL;
  }

  /*
   *      Remove all cell2 related arrays needed during collision handling
   */
  if (home->cell2) {
    free(home->cell2);
    home->cell2 = NULL;
  }
  home->cellCount = 0;
  home->nativeCellCount = 0;
  home->cell2nx = 0;
  home->cell2ny = 0;
  home->cell2nz = 0;

  if (home->cell2QentArray) {
    free(home->cell2QentArray);
    home->cell2QentArray = NULL;
  }

  /*
   *      Free the buffer used to hold the global cell charge tensor
   *      needed for calculating the far field stresses when the fast
   *      multipole code is not enabled.
   */
  if (home->cellCharge) {
    free(home->cellCharge);
    home->cellCharge = NULL;
  }

  /*
   *      Free all memory associated with the domain decomposition
   */
  if (home->decomp != (void *)NULL) {
    FreeDecomp(home, home->decomp);
    home->decomp = (void *)NULL;
  }

  /*
   *      Free arrays used in control and data file parsing
   */
  if (home->ctrlParamList != (ParamList_t *)NULL) {
    if (home->ctrlParamList->paramCnt > 0) {
      free(home->ctrlParamList->varList);
    }
    free(home->ctrlParamList);
  }
  home->ctrlParamList = NULL;

  if (home->dataParamList != (ParamList_t *)NULL) {
    if (home->dataParamList->paramCnt > 0) {
      free(home->dataParamList->varList);
    }
    free(home->dataParamList);
  }
  home->dataParamList = NULL;

  if (home->stressParamList != (ParamList_t *)NULL) {
    if (home->stressParamList->paramCnt > 0) {
      free(home->stressParamList->varList);
    }
    free(home->stressParamList);
  }
  home->stressParamList = NULL;

  /*
   *      Remove the burgers vector list and associated plane list
   *      if they were allocated.
   */
  if (home->burgData.burgList != (real8(*)[3])NULL) {
    free(home->burgData.burgList);
    home->burgData.burgList = (real8(*)[3])NULL;
  }

  if (home->burgData.planeList != (real8(*)[3])NULL) {
    free(home->burgData.planeList);
    home->burgData.planeList = (real8(*)[3])NULL;
  }

  if (home->burgData.numPlanesPerBurg != (int *)NULL) {
    free(home->burgData.numPlanesPerBurg);
    home->burgData.numPlanesPerBurg = (int *)NULL;
  }

  if (home->burgData.burgFirstPlaneIndex != (int *)NULL) {
    free(home->burgData.burgFirstPlaneIndex);
    home->burgData.burgFirstPlaneIndex = (int *)NULL;
  }

  free(home->param);
  home->param = NULL;
}

void ParadisFinish(Home_t *home) {
  int maxMem;

  if (home->myDomain == 0)
    printf("ParadisFinish\n");

  // Generate all type of output appropriate at program completion
  GenerateOutput(home, STAGE_TERM);

  /*
   *	Print out memory usage info for processor zero
   */
  if (home->myDomain == 0) {
    Meminfo(&maxMem);
    printf("Estimated memory usage (proc 0): %dk bytes\n", maxMem);
  }

  DSMPI::Finalize();

  ReleaseMemory(home);
  FreeCellCenters();
  if (home != NULL) {
    free(home);
  }
}
