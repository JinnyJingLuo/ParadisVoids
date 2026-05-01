/*****************************************************************************
 *
 *      Module:         Remesh.c
 *      Description:    This module contains functions common to more
 *                      than 1 of the supported version of remesh, plus
 *                      a generic entry function that invokes the
 *                      proper remesh version.
 *
 *      Included functions:
 *              CutSurfaceSegments()
 *              Remesh()
 *
 *****************************************************************************/
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include "Home.h"
#include "Comm.h"
#include "Topology.h"
#include "TwinPlaneCrossSlip.h"

#ifdef PARALLEL
#include "mpi.h"
#endif

/*-------------------------------------------------------------------------
 *
 *      Function:       Remesh
 *      Description:    This is just a generic function to set up
 *                      for a remesh operation then select and execute
 *                      the proper remesh function.
 *
 *-----------------------------------------------------------------------*/
void Remesh(Home_t *home) {
  int i;
  Node_t *node;
  Param_t *param;

  param = home->param;

#ifdef PARALLEL
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  ClearOpList(home);
  InitTopologyExemptions(home);
  RemeshRule_2(home);
  // junjie
  RemeshCSplane(home);
  // junjie

  /*
   *      Send to the neighboring domains, a list of all local
   *      operations (add/delete nodes, relinking of arms,
   *      etc) that may affect nodes in remote domains, and process
   *      the remesh operations from the neighboring domains
   */
  CommSendRemesh(home);
  FixRemesh(home);

  /*
   *      Under certain circumstances, parallel topological changes
   *      can create double links between nodes; links which can not be
   *      detected until after FixRemesh() is called... so, a quick
   *      check has to be done to clean up these potential double-
   *      links here, or they will cause problems later on.  Should
   *      only have to check nodes local to this domain.
   */
  for (i = 0; i < home->newNodeKeyPtr; i++) {
    if ((node = home->nodeKeys[i]) == (Node_t *)NULL)
      continue;
    (void)RemoveDoubleLinks(home, node, 0);
    node->flags &= ~NODE_CHK_DBL_LINK;
  }

  /*
   *      It is possible that remesh/coarsen has left an orphaned node.
   *      We need to get rid of any such nodes before the exchange of
   *      node information in CommSendGhosts().
   */
  RemoveOrphanedNodes(home);

#if PARALLEL
  MPI_Barrier(MPI_COMM_WORLD);
#endif
}
