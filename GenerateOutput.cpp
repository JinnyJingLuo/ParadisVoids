/*-------------------------------------------------------------------------
 *
 *      Module:      GenerateOutput.c
 *      Description:
 *
 *      Includes public functions:
 *
 *          DoParallelIO()
 *          GenerateOutput()
 *          GetOutputTypes()
 *          GetParIOGroup()
 *
 *------------------------------------------------------------------------*/
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "Home.h"
#include "Node.h"
#include "Comm.h"
#include "Util.h"
#include "WriteProp.h"
#include "Restart.h"

#ifdef PARALLEL
#include <mpi.h>
#endif

/*-------------------------------------------------------------------------
 *
 *      Function:    SendWriteToken
 *      Description: Send the write token to the specified task
 *
 *      Arguments:
 *          taskRank Rank of the task (in MPI_COMM_WORLD) to which
 *                   the write token is to be sent.
 *
 *------------------------------------------------------------------------*/
static void SendWriteToken(int taskRank) {
#ifdef PARALLEL
  int writeToken = 0;

  MPI_Send(&writeToken, 1, MPI_INT, taskRank, MSG_TOKEN_RING, MPI_COMM_WORLD);
#endif
  return;
}

/*-------------------------------------------------------------------------
 *
 *      Function:    RecvWriteToken
 *      Description: Wait for the write token from the specified task
 *
 *      Arguments:
 *          taskRank Rank of the task (in MPI_COMM_WORLD) from which
 *                   the write token is expected.
 *
 *------------------------------------------------------------------------*/
static void RecvWriteToken(int taskRank) {
#ifdef PARALLEL
  int writeToken = 0;
  MPI_Status reqStatus;

  MPI_Recv(&writeToken, 1, MPI_INT, taskRank, MSG_TOKEN_RING, MPI_COMM_WORLD,
           &reqStatus);
#endif
  return;
}

/*-------------------------------------------------------------------------
 *
 *      Function:    GetOutputTypes
 *      Description: Given the stage of program execution, determine
 *                   which types of output need to be generated and
 *                   increment any applicable file sequence counters for
 *                   those output types. -- Note the counters are not to
 *                   be updated while dumping data during program termination.
 *      Args:
 *          stage        indicates the program execution stage which
 *                       determines the types of output this function
 *                       might produce.  Valid values for this field are:
 *
 *                               STAGE_INIT
 *                               STAGE_CYCLE
 *                               STAGE_TERM
 *
 *          outputTypes  pointer to an integer which on exit from this
 *                       subroutine will contain a bitmask with 1 bit
 *                       set for each type of output to be generated.
 *
 *------------------------------------------------------------------------*/
static void GetOutputTypes(Home_t *home, int stage, int *outputTypes) {
  real8 timeNow, dumpTime;
  Param_t *param;

  static int xWinAlive = 1;

  param = home->param;
  timeNow = param->timeNow;

  /*
   *      If the code is
   *      in the termination stage, or it is at the end of a cycle
   *      and either the elapsed simulation time since the last dump
   *      has exceeded the allowable delta time between restart dumps
   *      OR the current cycle is a multiple of the specified restart
   *      dump frequency.
   *
   *      NOTE: The last two conditions are mutually exclusive. If
   *            a delta time has been specified to determine the
   *            time to dump files, no check of the current
   *            cycle against the dump frequency will be done.
   */
  if (stage == STAGE_TERM) {
    *outputTypes |= GEN_RESTART_DATA;
  } else if ((stage == STAGE_CYCLE) &&
             ((home->cycle % param->savefreq) ==
              0)) // qjiao: for debug to controll the exact step for output.
  // else if ((stage == STAGE_CYCLE) && ((home->cycle >= 16000) && (home->cycle
  // <= 18000)))
  {
    *outputTypes |= GEN_RESTART_DATA;
  }

  /*
   *      If TECPLOT file dumps are enabled, do so if the code is
   *      in the termination stage, or it is at the end of a cycle
   *      and either the elapsed simulation time since the last dump
   *      has exceeded the allowable delta time between file dumps
   *      OR the current cycle is a multiple of the specified
   *      dump frequency.
   *
   *      NOTE: The last two conditions are mutually exclusive. If
   *            a delta time has been specified to determine the
   *            time to dump files, no check of the current
   *            cycle against the dump frequency will be done.
   */
  if (param->tecplot) {
    if (stage == STAGE_TERM) {
      *outputTypes |= GEN_TECPLOT_DATA;
    } else if (stage == STAGE_CYCLE) {
      if (param->tecplotdt > 0.0) {
        dumpTime = param->tecplottime + param->tecplotdt;
        if (timeNow >= dumpTime) {
          param->tecplottime = timeNow;
          param->tecplotcounter++;
          *outputTypes |= GEN_TECPLOT_DATA;
        }
      } else if ((home->cycle % param->savefreq) == 0) {
        param->tecplotcounter++;
        *outputTypes |= GEN_TECPLOT_DATA;
      }
    }
  }

  /*
   *      Check if creation of the density field file is enabled.  Only
   *      done at code termination time...
   */
  if ((stage == STAGE_TERM) &&
      ((param->savedensityspec[0] > 0) && (param->savedensityspec[1] > 0) &&
       (param->savedensityspec[2] > 0))) {
    *outputTypes |= GEN_DENSITY_DATA;
  }

  /*
   *      If properties file dumps are enabled, do so if the code is
   *      at the end of a cycle and either the elapsed simulation
   *      time since the last dump has exceeded the allowable delta
   *      time between file dumps OR the current cycle is a multiple
   *      of the specified dump frequency.
   *
   *      NOTE: The last two conditions are mutually exclusive. If
   *            a delta time has been specified to determine the
   *            time to dump files, no check of the current
   *            cycle against the dump frequency will be done.
   */
  if (stage == STAGE_CYCLE) {
    if (param->savepropdt > 0.0) {
      dumpTime = param->saveproptime + param->savepropdt;
      if (timeNow >= dumpTime) {
        param->saveproptime = timeNow;
        *outputTypes |= GEN_PROPERTIES_DATA;
      }
    } else if ((home->cycle % param->savefreq) == 0) {
      *outputTypes |= GEN_PROPERTIES_DATA;
    }
  }
}

/*-------------------------------------------------------------------------
 *
 *      Function:    GetParIOGroup
 *      Description: Determine the parallel I/O group of which the current
 *                   task is a member, as well as the previous and
 *                   next tasks in the same group, etc.
 *
 *------------------------------------------------------------------------*/
void GetParallelIOGroup(Home_t *home) {
  int numTasks, thisTask;
  int taskID, group;
  int numIOGroups;
  int smallGroupSize, largeGroupSize, thisGroupSize;
  int numLargeGroups;
  int firstInGroup, lastInGroup, thisGroup;
  int *taskList;
  Param_t *param;
#ifdef PARALLEL
  MPI_Group groupCommWorld, groupLastInIOGroup;
#endif

  param = home->param;

  numIOGroups = param->numIOGroups;
  numTasks = home->numDomains;
  thisTask = home->myDomain;

  smallGroupSize = numTasks / numIOGroups;
  numLargeGroups = numTasks % numIOGroups;
  largeGroupSize = smallGroupSize + (numLargeGroups > 0);

  taskList = (int *)malloc(numIOGroups * sizeof(int));

  /*
   *      First make a list including only the processes that are
   *      the last in their respective IO groups
   */

  for (taskID = 0, group = 0; group < numIOGroups; group++) {
    thisGroupSize = (group < numLargeGroups) ? largeGroupSize : smallGroupSize;

    lastInGroup = taskID + (thisGroupSize - 1);

    taskList[group] = lastInGroup;
    taskID += thisGroupSize;
  }

  /*
   *      Now calculate the task specific stuff like group ID, first
   *      and last tasks in the group, and so on.
   */
  if ((numLargeGroups * largeGroupSize) > thisTask) {
    thisGroup = thisTask / largeGroupSize;
    firstInGroup = thisGroup * largeGroupSize;
    lastInGroup = firstInGroup + largeGroupSize - 1;
  } else {
    thisGroup =
        numLargeGroups +
        ((thisTask - (numLargeGroups * largeGroupSize)) / smallGroupSize);
    firstInGroup = numLargeGroups * largeGroupSize +
                   (thisGroup - numLargeGroups) * smallGroupSize;
    lastInGroup = firstInGroup + smallGroupSize - 1;
  }

  home->ioGroupNum = thisGroup;
  home->firstInIOGroup = firstInGroup;
  home->lastInIOGroup = lastInGroup;
  home->prevInIOGroup = MAX(thisTask - 1, firstInGroup);
  home->nextInIOGroup = MIN(thisTask + 1, lastInGroup);
  home->isLastInIOGroup = (thisTask == lastInGroup);
  home->isFirstInIOGroup = (thisTask == firstInGroup);

#ifdef PARALLEL
  /*
   *      Create a new MPI group containing the processes identified above,
   *      and a new communicator encompasing those processes.
   */

  MPI_Comm_group(MPI_COMM_WORLD, &groupCommWorld);
  MPI_Group_incl(groupCommWorld, numIOGroups, taskList, &groupLastInIOGroup);
  MPI_Comm_create(MPI_COMM_WORLD, groupLastInIOGroup, &home->commLastInIOGroup);

  /*
   *      Free the MPI group definitions which are no longer needed
   */
  MPI_Group_free(&groupCommWorld);
  MPI_Group_free(&groupLastInIOGroup);

#endif

  free(taskList);
}

/*-------------------------------------------------------------------------
 *
 *      Function:    GetCounts
 *
 *      Description: Obtain both local and global node/segment counts.
 *
 *      Args:
 *          numNodes      location in which to return to caller the total
 *                        number of nodes in the simulation.
 *          numSegs       location in which to return to caller the total
 *                        number of unique segments in the simulation.
 *          numArms       location in which to return the sum arm count
 *                        of all nodes in the simulation
 *          numLocalNodes location in which to return the number of nodes
 *                        in the local domain.
 *          numLocalSeg   location in which to return the number of unique
 *                        segments in the local domain.
 *
 *      Last Modified: 09/10/2008 - gh  Modified to add the <numLocalNodes>
 *                                  and <numLocalSegs> parameters.
 *
 *------------------------------------------------------------------------*/
static void GetCounts(Home_t *home, int *numNodes, int *numSegs, int *numArms,
                      int *numLocalNodes, int *numLocalSegs) {
  int i, j;
  int locVals[3], globalVals[3];
  Node_t *node;

  locVals[0] = 0;
  locVals[1] = 0;
  locVals[2] = 0;

  globalVals[0] = 0;
  globalVals[1] = 0;
  globalVals[2] = 0;

  for (i = 0; i < home->newNodeKeyPtr; i++) {

    if ((node = home->nodeKeys[i]) == (Node_t *)NULL) {
      continue;
    }

    locVals[0] += 1;
    locVals[2] += node->numNbrs;

    for (j = 0; j < node->numNbrs; j++) {
      if ((home->myDomain == node->nbrTag[j].domainID) &&
          (node->nbrTag[j].index < i)) {
        continue;
      }
      locVals[1] += 1;
    }
  }

#ifdef PARALLEL
  MPI_Reduce(locVals, globalVals, 3, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

  *numNodes = globalVals[0];
  *numSegs = globalVals[1];
  *numArms = globalVals[2];
#else
  *numNodes = locVals[0];
  *numSegs = locVals[1];
  *numArms = locVals[2];
#endif

  *numLocalNodes = locVals[0];
  *numLocalSegs = locVals[1];

  return;
}

/*-------------------------------------------------------------------------
 *
 *      Function:    DoParallelIO
 *      Description:
 *
 *      Args:
 *          outputTypes  integer bitfield.  All non-zero bits correspond
 *                       to specific types of output that are to be
 *                       generated at this stage of program execution.
 *                       See Util.h for the mapping of output types
 *                       to bit positions.
 *          stage        indicates the program execution stage which
 *                       determines the types of output this function
 *                       might produce.
 *
 *------------------------------------------------------------------------*/
static void DoParallelIO(Home_t *home, int outputTypes, int stage) {
  int ioGroup, prevInGroup, nextInGroup, numIOGroups;
  int thisDomain, isFirstInGroup, isLastInGroup;
  int writePrologue, writeEpilogue;
  int sendToken, recvToken;
  int writeToken = 0, numSegs = 0, numArms = 0, totFragmentCount = 0;
  int nodesWritten = 0, segsWritten = 0;
  int countInGroup[2] = {0, 0};
  char baseName[128];
  void *fragmentList;
  Param_t *param;
  BinFileData_t binData;
#ifdef PARALLEL
  MPI_Request req;
  MPI_Status reqStatus;
#endif

  thisDomain = home->myDomain;
  param = home->param;

  numIOGroups = param->numIOGroups;

  memset(&binData, 0, sizeof(BinFileData_t));

  binData.firstInGroup = home->firstInIOGroup;
  binData.lastInGroup = home->lastInIOGroup;

  ioGroup = home->ioGroupNum;

  isFirstInGroup = home->isFirstInIOGroup;
  isLastInGroup = home->isLastInIOGroup;

  prevInGroup = home->prevInIOGroup;
  nextInGroup = home->nextInIOGroup;

  recvToken = (thisDomain != prevInGroup);
  sendToken = (thisDomain != nextInGroup);

  /*
   *      Domain zero (first member of first I/O group) will always
   *      do output file creation and initialization, and the last
   *      member of the last I/O group always adds any necessary
   *      trailers/whatever to the file
   */
  writePrologue = (thisDomain == 0);
  writeEpilogue = ((ioGroup == (numIOGroups - 1)) && isLastInGroup);

  /*
   *      Certain output routines require total counts of nodes,
   *      segments, etc.  If we're doing any of these types of output
   *      then do a global operation to get these values before
   *      beginning the I/O.
   */
  if ((outputTypes & GEN_RESTART_DATA) || (outputTypes & GEN_TECPLOT_DATA)) {
    GetCounts(home, &param->nodeCount, &numSegs, &numArms, &binData.nodeCount,
              &binData.segCount);
  }

  /*
   *      Call all I/O functions...
   *
   *      Before writing any specific output file, wait for the write
   *      token from the predecesor in the I/O group.  When done
   *      write the specific output type, if this task is not the
   *      last in its I/O group, then send the write token to the
   *      next process in the I/O group so the other task can start
   *      on that output type while the current task continues on
   *      with other output types.
   */
  if ((outputTypes & GEN_RESTART_DATA) != 0) {
    if (stage == STAGE_CYCLE) {
      snprintf(baseName, sizeof(baseName), "rs%04d", param->savecounter);
    } else {
      snprintf(baseName, sizeof(baseName), "restart.cn");
    }
    if (recvToken)
      RecvWriteToken(prevInGroup);
    WriteRestart(home, baseName, ioGroup, isFirstInGroup, writePrologue,
                 writeEpilogue);
    if (sendToken)
      SendWriteToken(nextInGroup);
  }

  if ((outputTypes & GEN_TECPLOT_DATA) != 0) {
    if (stage == STAGE_CYCLE) {
      snprintf(baseName, sizeof(baseName), "tecdata%04d",
               param->tecplotcounter);
    } else {
      snprintf(baseName, sizeof(baseName), "tecdata.final");
    }
    if (recvToken)
      RecvWriteToken(prevInGroup);
    Tecplot(home, baseName, ioGroup, isFirstInGroup, writePrologue,
            writeEpilogue, numSegs);
    if (sendToken)
      SendWriteToken(nextInGroup);
  }

  /*
   *      If we wrote the text restart file we still need to store the
   *      name of the recently written restart file to disk.  This
   *      involves an explicit syncronization point in the code since
   *      we don't want to do this until all processes have completed
   *      writing their restart data.
   */
  if ((outputTypes & GEN_RESTART_DATA) != 0) {
    if (stage == STAGE_CYCLE) {
      snprintf(baseName, sizeof(baseName), "rs%04d", param->savecounter);
    } else {
      snprintf(baseName, sizeof(baseName), "restart.cn");
    }
#ifdef PARALLEL
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    if (thisDomain == 0) {
      SetLatestRestart(baseName);
    }
  }
}

/*-------------------------------------------------------------------------
 *
 *      Function:    GenerateOutput
 *      Description: This subroutine controls generation of all types
 *                   of output appropriate to the current stage of program
 *                   execution.
 *      Args:
 *          stage    indicates the program execution stage which determines
 *                   the types of output this function might produce.  Valid
 *                   values for this field are:
 *
 *                           STAGE_INIT
 *                           STAGE_CYCLE
 *                           STAGE_TERM
 *
 *------------------------------------------------------------------------*/
void GenerateOutput(Home_t *home, int stage) {
  int outputTypes = 0;
  real8 localDensity, globalDensity;
  char fileName[128];
  time_t tp;
  Param_t *param;

  param = home->param;

  /*
   *      Determine what types of output need to be generated at
   *      this stage.
   */
  GetOutputTypes(home, stage, &outputTypes);

  /*
   *      Dump slip amount and strain decomposition?
   */

  if (outputTypes & GEN_DENSITY_DATA) {
    WriteDensityField(home, "densityfield.out");
  }
  /*
   *      Preserve various properites i.e. density, etc.
   *               (property vs. time files)
   */
  if ((outputTypes & GEN_PROPERTIES_DATA) != 0) {

    WriteProp(home, DENSITY); /* Actually, length of disloctions   */
                              /* both total and broken down by     */
                              /* groups of burgers vector types    */

    WriteProp(home, DENSITY_DELTA); /* Per-burgers vector density */
                                    /* gain/loss                  */

    WriteProp(home, EPSDOT); /* Resultant strain rate (scalar)    */

    WriteProp(home, ALL_EPS); /* Plastic strain tensor             */

    WriteProp(home, EPS); /* Resultant plastic strain (scalar) */
  }

  /*
   *      A number of the types of output will be done (potentially)
   *      in parallel.  Handle those in a separate routine
   */
  if (((outputTypes & GEN_RESTART_DATA) != 0) ||
      ((outputTypes & GEN_TECPLOT_DATA) != 0)) {
    DoParallelIO(home, outputTypes, stage);
  }

#ifdef PARALLEL
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  // finally, if there are shearable precipitates, write the APB data
  if (home->param->ShearablePrecipitates != 0) {
    if (home->cycle % home->param->savefreq == 0) {
      home->poPrecipitateServer->WriteAPB(string(home->param->APBFileName));
    }
  }
}
