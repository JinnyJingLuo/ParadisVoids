/**************************************************************************
 *
 *      Module:       ReadRestart.c
 *      Description:  This module contains the functions for reading
 *                    parameters and nodal data from the control and
 *                    nodal data files indicated (or implied) by the
 *                    command line arguments.
 *
 *                    Several versions of the nodal data file are
 *                    supported, so there are multiple functions
 *                    provided to deal with different formats.
 *
 *                    NOTE:  The paradisconvert utility can be used
 *                    to translate any old style control and nodal
 *                    data files to the most current version.
 *
 *      Included public functions:
 *          AssignNodesToDomains()
 *          FreeNodeLists()
 *          ReadControlFile()
 *          ReadNodeDataFile()
 *
 *      Included private functions:
 *          ReadPreV4DataParams()
 *
 *************************************************************************/
#include "Home.h"
#include "InData.h"
#include "Tag.h"
#include "Util.h"
#include "Decomp.h"
#include "Restart.h"
#include "Parse.h"
#include "DSMPI.h"

/*---------------------------------------------------------------------------
 *
 *      Function:       ReadControlFile
 *      Description:    Read and parse the contents of the control
 *                      parameter file, saving the values associated
 *                      with the parameters in the corresponding
 *                      variables.
 *
 *      Arguments:
 *          ctrlFileName  Name of the control file to be read.
 *
 *-------------------------------------------------------------------------*/
void ReadControlFile(Home_t *home, const char *ctrlFileName) {
  int maxTokenLen, tokenType, pIndex, okay;
  int valType, numVals;
  char token[256];
  void *valList;
  FILE *fpCtrl;

  maxTokenLen = sizeof(token);
  tokenType = TOKEN_GENERIC;

  if ((fpCtrl = fopen(ctrlFileName, "r")) == (FILE *)NULL) {
    Fatal("ReadControlFile: Error %d opening file %s", errno, ctrlFileName);
  }

  while (1) {

    /*
     *          Next token (if any) should be a parameter name...
     */
    tokenType = GetNextToken(fpCtrl, token, maxTokenLen);

    /*
     *          If we hit the end of the file (i.e. no more tokens)
     *          we're done...
     */
    if (tokenType == TOKEN_NULL) {
      break;
    }

    /*
     *          If there were any parsing errors, just abort
     */
    if (tokenType == TOKEN_ERR) {
      Fatal("ReadControlFile: Error parsing file %s", ctrlFileName);
    }

    /*
     *          Obtain values associated with the parameter; don't
     *          save values for unknown/obsolete parameters.
     */
    pIndex = LookupParam(home->ctrlParamList, token);

    if (pIndex < 0) {
      valType = V_NULL;
      numVals = 0;
      valList = (void *)NULL;
    } else {
      valType = home->ctrlParamList->varList[pIndex].valType;
      numVals = home->ctrlParamList->varList[pIndex].valCnt;
      valList = home->ctrlParamList->varList[pIndex].valList;
      home->ctrlParamList->varList[pIndex].flags |= VFLAG_SET_BY_USER;
    }

    okay = GetParamVals(fpCtrl, valType, numVals, valList);
    if (!okay) {
      Fatal("Error obtaining values for parameter %s",
            home->ctrlParamList->varList[pIndex].varName);
    }
  }

  fclose(fpCtrl);

  return;
}

/*---------------------------------------------------------------------------
 *
 *      Function:       AssignNodesToDomains
 *      Description:    Using the provided domain decomposition data,
 *                      loop though all nodes in the inData node array
 *                      and assign each node to the domain which encompasses
 *                      the node's coordinates.
 *
 *      Arguments:
 *          nodeCount   number of nodes in the inData node array.
 *          nodeLists   Location in which to return to the caller the
 *                      array of lists of nodes to be sent to the
 *                      remote domains.  The nodeLists array has
 *                      a list pointer for every domain in the problem.
 *                      (each list is an array of indices into the
 *                      inData node array)
 *          listCounts  Location in which to return to the caller an
 *                      array of integers indicating the count of nodes
 *                      on the corresponding list in <nodeLists>.
 *
 *-------------------------------------------------------------------------*/
void AssignNodesToDomains(Home_t *home, InData_t *inData, int ***nodeLists,
                          int **listCounts) {
  int i, j;
  int nXdoms, nYdoms, nZdoms;
  int domIndex, nDoms, len, nexti;
  int *qHead, *qCount, *list;
  int **listArray, *countArray;
  Param_t *param;
  Node_t *node;

  param = home->param;
  int nodeCount = param->nodeCount;
  nXdoms = param->nXdoms;
  nYdoms = param->nYdoms;
  nZdoms = param->nZdoms;

  if (nodeCount == 0) {
    *nodeLists = (int **)NULL;
    *listCounts = (int *)NULL;
    return;
  }

  nDoms = nXdoms * nYdoms * nZdoms;

  listArray = (int **)calloc(1, nDoms * sizeof(int *));
  countArray = (int *)calloc(1, nDoms * sizeof(int));

  /*
   *      Allocate and initialize an array of integers (1 per domain)
   *      to be used as pointers to the head of the queue of nodes
   *      assigned to the domains.
   */
  qHead = (int *)malloc(nDoms * sizeof(int));
  qCount = (int *)malloc(nDoms * sizeof(int));

  for (i = 0; i < nDoms; i++) {
    qHead[i] = -1;
    qCount[i] = 0;
  }

  /*
   *      Loop through all the nodes on the current node list, find
   *      the proper domain for the node based on the node coordinates
   *      and add the node to the domain's queue.
   */
  for (i = 0; i < nodeCount; i++) {

    node = &inData->node[i];
    domIndex = FindCoordDomain(home, node->x, node->y, node->z);

    /*
     *          Each node on a queue contains the index of the next node
     *          on the queue.  To add a node to the queue, we just set
     *          that node's 'next node' pointer to the current head of
     *          the queue, and then point the queue head to the new node.
     */
    nexti = qHead[domIndex];
    if (nexti < 0) {
      node->next = (Node_t *)NULL;
    } else {
      node->next = &inData->node[qHead[domIndex]];
    }
    qHead[domIndex] = i;
    qCount[domIndex]++;
  }

  /*
   *      For each domain, generate the list of indices in the node array
   *      for nodes to be sent to the domain.
   */
  for (domIndex = 0; domIndex < nDoms; domIndex++) {

    list = (int *)NULL;
    len = qCount[domIndex];
    nexti = qHead[domIndex];

    if (len > 0)
      list = (int *)malloc(len * sizeof(int));

    for (j = 0; j < len; j++) {
      list[j] = nexti;
      /*
       *              Do some pointer arithmetic to get convert pointer
       *              addresses into array indices.
       */
      if (inData->node[nexti].next == (Node_t *)NULL) {
        nexti = -1;
      } else {
        nexti = inData->node[nexti].next - inData->node;
      }
    }

    if (nexti != -1) {
      Fatal("Queue error. domain %d, queue len %d", domIndex, len);
    }

    listArray[domIndex] = list;
    countArray[domIndex] = len;
  }

  /*
   *      Free up the now unneeded arrays.  (the lists of node indices
   *      will be freed elsewhere when they are no longer needed)
   */
  free(qHead);
  free(qCount);

  *nodeLists = listArray;
  *listCounts = countArray;
}

/*---------------------------------------------------------------------------
 *
 *      Function:       FreeNodeLists
 *      Description:    Releases temporary storage associated with
 *                      the arrays of node indices identifying nodes
 *                      to be sent to remote domains.
 *
 *      Arguments:
 *          nodeLists   Address of the the array of lists of nodes to
 *                      that were sent to the remote domains.  The nodeLists
 *                      array has a list pointer for every domain in the
 *                      problem.  On return, the contents of this adress
 *                      will be zeroed.
 *          listCounts  Address of the array of integers indicating the
 *                      count of nodes on the corresponding list in
 *                      <nodeLists>.  On return, the contents of this
 *                      adress will be zeroed.
 *
 *-------------------------------------------------------------------------*/
void FreeNodeLists(Home_t *home, int ***nodeLists, int **listCounts) {
  int dom;

  if (*listCounts != (int *)NULL) {
    free(*listCounts);
    *listCounts = (int *)NULL;
  }

  if (*nodeLists != (int **)NULL) {
    for (dom = 0; dom < home->numDomains; dom++) {
      if ((*nodeLists)[dom] != (int *)NULL) {
        free((*nodeLists)[dom]);
        (*nodeLists)[dom] = (int *)NULL;
      }
    }
    free(*nodeLists);
    *nodeLists = (int **)NULL;
  }
}

/*------------------------------------------------------------------------
 *
 *	Function:	ReadNodeDataFile
 *	Description:	Read the nodal data from a single or multi-segment
 *                      data file, assign the nodes to appropriate domains
 *                      based on the node coordinates and the current domain
 *			decomposition, and distribute the nodal data
 *			to the appropriate remote domains.
 *
 *			The data will be processed in blocks so as to
 *			avoid loading the entire problem space onto
 *			processor zero and potentially exhausting the
 *			memory on small-memory systems.
 *
 *			NOTE: This function is only for use with newer
 *			data files, not for reading node data provided
 *			in the old format within the control file itself!
 *
 *      Arguments:
 *          inData    pointer to structure in which to temporarily
 *                    store domain decomposition, nodal data, etc
 *                    from the restart file.
 *          dataFile  Name of the nodal data file.  For segmented
 *                    restart files, this may be the base name
 *                    (i.e. no trailing sequence number) or the
 *                    name of the first file segment.
 *
 *-----------------------------------------------------------------------*/
void ReadNodeDataFile(Home_t *home, InData_t *inData, const char *dataFile) {
  int dom, iNbr;
  int maxTokenLen, tokenType, pIndex, okay;
  int valType, numVals;
  int numDomains;
  int numNbrs;
  int fileSeqNum = 0;
  int *globalMsgCnt, *localMsgCnt, **nodeLists, *listCounts;
  int localNodeCount, globalNodeCount;
  int binRead = 0;
  int nextAvailableTag = 0;
  real8 burgSumX, burgSumY, burgSumZ;
  void *valList;
  char inLine[500], token[256];
  char baseFileName[256], tmpFileName[256];
  FILE *fpSeg;
  Node_t *node;
  Param_t *param;
  ParamList_t *dataParamList;

  param = home->param;
  numDomains = home->numDomains;
  fpSeg = (FILE *)NULL;

  dataParamList = home->dataParamList;
  globalMsgCnt = (int *)malloc((numDomains + 1) * sizeof(int));
  localMsgCnt = (int *)malloc((numDomains + 1) * sizeof(int));

  memset(inLine, 0, sizeof(inLine));
  maxTokenLen = sizeof(token);

  /*
   *      Only domain zero reads the initial stuff...
   */
  if (home->myDomain == 0) {
    if (dataFile == NULL) {
      Fatal("ReadNodeDataFile(): ReadNodeDataFile: No data file provided");
    }
    snprintf(tmpFileName, sizeof(tmpFileName), "%s", dataFile);
    if (!strcmp(&tmpFileName[strlen(tmpFileName) - 2], ".0")) {
      tmpFileName[strlen(tmpFileName) - 2] = 0;
    }
    snprintf(baseFileName, sizeof(baseFileName), "%s", tmpFileName);

    /*
     *          Try to open the nodal data file.  If the specified file
     *          can not be opened, try looking for a segmented data file
     *          (i.e. look for the first file segment <dataFile>.0). If
     *          that can't be opened either, then exit with an error.
     */

    if ((fpSeg = fopen(tmpFileName, "r")) == NULL) {
      snprintf(tmpFileName, sizeof(tmpFileName), "%s.%d", baseFileName,
               fileSeqNum);
      if ((fpSeg = fopen(tmpFileName, "r")) == NULL) {
        Fatal("ReadNodeDataFile(): Error %d opening file %s to read nodal data",
              errno, dataFile);
      }
    }

    /*
     *          Get the first token.  This should either be a known
     *          parameter identifier, or a file version number.
     */
    tokenType = GetNextToken(fpSeg, token, maxTokenLen);
    pIndex = LookupParam(dataParamList, token);

    /*
     *              Just go through the nodal data file reading all
     *              the associated parameters.  Need to do special
     *              processing of the domain decomposition, and when
     *              when we hit the 'nodaldata' identifier, just break
     *              out of this loop.
     */
    while ((tokenType != TOKEN_ERR) && (tokenType != TOKEN_NULL)) {
      if (pIndex >= 0) {
        /*
         *                      Token represents a known parameter identifier,
         * so read the associated value(s).
         */
        valType = dataParamList->varList[pIndex].valType;
        numVals = dataParamList->varList[pIndex].valCnt;
        valList = dataParamList->varList[pIndex].valList;
        okay = GetParamVals(fpSeg, valType, numVals, valList);
        if (!okay) {
          Fatal("ReadNodeDataFile(): Parsing Error obtaining values for "
                "parameter %s\n",
                dataParamList->varList[pIndex].varName);
        }
      } else {
        /*
         *                      Token does not represent one of the simple
         *                      parameters.  If it's not one of the identifiers
         *                      that needs special handling, skip it.
         */
        if (strcmp(token, "nodalData") == 0) {
          /*
           *                          When we hit the nodal data, we can break
           *                          out of the loop because we are assuming
           *                          all other data file parameters have been
           *                          processed.  If they have not, we have a
           *                          problem since processing of the nodal data
           *                          requires the other parameters.
           *                          Note: Remainder of the file should just
           *                          contain " = " followed by the nodal data,
           *                          so be sure to skip the next token before
           *                          processing the nodal data.
           */
          tokenType = GetNextToken(fpSeg, token, maxTokenLen);
          break;
        }
      }
      tokenType = GetNextToken(fpSeg, token, maxTokenLen);
      if ((tokenType == TOKEN_NULL) || (tokenType == TOKEN_ERR)) {
        Fatal("ReadNodeDataFile(): Parsing error on file %s\n", tmpFileName);
      }
      pIndex = LookupParam(dataParamList, token);
    }

    /*
     *          Need to set some of the values that are dependent on
     *          the simulation size before we go any further.
     */

    param->Lx = param->Dimensions[X];
    param->Ly = param->Dimensions[Y];
    param->Lz = param->Dimensions[Z];

    param->invLx = 1.0 / param->Lx;
    param->invLy = 1.0 / param->Ly;
    param->invLz = 1.0 / param->Lz;
  } /* if (home->myDomain == 0) */

  /*
   *      Domain zero now needs to pass various pieces of data to
   *      the remote domains before any processes begin to read
   *      and distribute the nodal data... start with the param structure.
   */
  DSMPI::Broadcast((char *)param, sizeof(Param_t), MPI_CHAR, 0);

  /*
   *      Now the base nodal data file name.
   */
  DSMPI::Broadcast((char *)baseFileName, sizeof(baseFileName), MPI_CHAR, 0);

  // some checks on the axes
  double dTolerance = 1.0E-6;
  Vector oXAxis(home->param->XAxis[0], home->param->XAxis[1],
                home->param->XAxis[2]);
  if (oXAxis.Length() < dTolerance) {
    Fatal("ReadNodeDataFile(): zero length X axis\n");
  }
  Vector oYAxis(home->param->YAxis[0], home->param->YAxis[1],
                home->param->YAxis[2]);
  if (oYAxis.Length() < dTolerance) {
    Fatal("ReadNodeDataFile(): zero length Y axis\n");
  }
  if (fabs(oXAxis * oYAxis) > dTolerance) {
    Fatal("ReadNodeDataFile(): X and Y axes are not orthogonal\n");
  }
  oXAxis.Normalize();
  oYAxis.Normalize();

  home->param->XAxis[0] = oXAxis.GetX();
  home->param->XAxis[1] = oXAxis.GetY();
  home->param->XAxis[2] = oXAxis.GetZ();

  home->param->YAxis[0] = oYAxis.GetX();
  home->param->YAxis[1] = oYAxis.GetY();
  home->param->YAxis[2] = oYAxis.GetZ();

  Vector oZAxis = oXAxis ^ oYAxis;
  if (oZAxis.Length() < dTolerance) {
    Fatal("ReadNodeDataFile(): zero length Z axis\n");
  }
  oZAxis.Normalize();

  home->param->ZAxis[0] = oZAxis.GetX();
  home->param->ZAxis[1] = oZAxis.GetY();
  home->param->ZAxis[2] = oZAxis.GetZ();

  // set the param of the input data object
  inData->param = home->param;
  /*
   *      All processes loop until all nodal data has been read in
   *      and distributed to the appropriate domains.
   */
  unsigned int i = 0;
  listCounts = NULL;
  nodeLists = NULL;

  // only domain 0 reads in the node data initially
  if (home->myDomain == 0) {
    if (fpSeg == NULL) {
      Fatal("ReadNodeDataFile(): NULL node file\n");
    }
    inData->node = (Node_t *)calloc(1, param->nodeCount * sizeof(Node_t));
    for (i = 0; i < param->nodeCount; i++) {
      Getline(inLine, sizeof(inLine), fpSeg);
      if ((inLine[0] == 0) && (i < param->nodeCount - 1)) {
        Fatal("ReadNodeDataFile(): node file ended before reading all nodes\n");
      }

      // get a new node
      node = &inData->node[i];
      // read node data
      if (home->param->EnableTwinPlaneCrossSlip == 1) {
        sscanf(inLine, "%d,%d %lf %lf %lf %d %d %d", &node->myTag.domainID,
               &node->myTag.index, &node->x, &node->y, &node->z, &numNbrs,
               &node->constraint, &node->csState);
      } else
        sscanf(inLine, "%d,%d %lf %lf %lf %d %d", &node->myTag.domainID,
               &node->myTag.index, &node->x, &node->y, &node->z, &numNbrs,
               &node->constraint);

      // make sure it is the primary image
      GetPrimaryImage(home->param, &node->x, &node->y, &node->z);
      node->oldx = node->x;
      node->oldy = node->y;
      node->oldz = node->z;
      node->VelocityDapmingSteps = 0;
      // read the surface normal
      if (home->param->dataFileVersion >= 5) {
        Getline(inLine, sizeof(inLine), fpSeg);
        sscanf(inLine, "\t\t%lf %lf %lf", &node->dNx, &node->dNy, &node->dNz);
      }
      AllocNodeArms(node, numNbrs);

      // read arm data
      burgSumX = 0.0;
      burgSumY = 0.0;
      burgSumZ = 0.0;

      for (iNbr = 0; iNbr < node->numNbrs; iNbr++) {
        Getline(inLine, sizeof(inLine), fpSeg);
        sscanf(inLine, "%d,%d %lf %lf %lf", &node->nbrTag[iNbr].domainID,
               &node->nbrTag[iNbr].index, &node->burgX[iNbr],
               &node->burgY[iNbr], &node->burgZ[iNbr]);
        Getline(inLine, sizeof(inLine), fpSeg);

        if (home->param->dataFileVersion >= 6) {
          sscanf(inLine, "%lf %lf %lf %d", &node->nx[iNbr], &node->ny[iNbr],
                 &node->nz[iNbr], &node->piChainID[iNbr]);
        } else {
          sscanf(inLine, "%lf %lf %lf", &node->nx[iNbr], &node->ny[iNbr],
                 &node->nz[iNbr]);
        }
        Normalize(&node->nx[iNbr], &node->ny[iNbr], &node->nz[iNbr]);

        burgSumX += node->burgX[iNbr];
        burgSumY += node->burgY[iNbr];
        burgSumZ += node->burgZ[iNbr];
      }

      /*
       *                  Just a quick sanity check, to make sure burgers
       *                  vector is conserved for all unconstrained nodes.
       */
      if (node->constraint == UNCONSTRAINED) {
        if ((fabs(burgSumX) > dTolerance) || (fabs(burgSumY) > dTolerance) ||
            (fabs(burgSumZ) > dTolerance)) {
          printf("Error: node (%d,%d)\n", node->myTag.domainID,
                 node->myTag.index);
          for (iNbr = 0; iNbr < node->numNbrs; iNbr++) {
            printf("  arm[%d] burg = %e %e %e\n", iNbr, node->burgX[iNbr],
                   node->burgY[iNbr], node->burgZ[iNbr]);
          }
          cout << "Coardiantes of node" << node->x << "\t" << node->y << "\t"
               << node->z << endl;
          Fatal("Burger's vector not conserved!");

          // cout <<"Coardiantes of neigbor" << nbr->x <<"\t" <<nbr->y <<"\t"
          // <<nbr->z<<endl; cout <<"IDs" <<"\t" << nbrArmID << "\t" <<
          // nbr->numNbrs <<endl;
        }
      }
    }
    // now since we have read all the node data, update the domain decomposition
    // so as to balance the number of nodes over the domains
    printf("Generating node-based decomposition ... ");
    inData->decomp = PointBasedDecomp(home, inData);
    printf("done\n");
  }
  BroadcastDecomp(home, inData->decomp);
  // distribute the nodes over the domains
  if (home->myDomain == 0) {
    AssignNodesToDomains(home, inData, &nodeLists, &listCounts);
  }

  /*
   *          Set up an array (1 element per domain).  Each reader
   *          task sets to 1 the array entry for each remote domain
   *          to which it will be sending data.  When we do a global
   *          reduction to sum up the arrays, the resulting array
   *          contains the number of messages that will be sent to
   *          each domain during this communication.
   *
   *          Plus 1 extra element set to 1 if ANY process is sending
   */
  memset(globalMsgCnt, 0, (numDomains + 1) * sizeof(int));
  memset(localMsgCnt, 0, (numDomains + 1) * sizeof(int));

  // all the domains will skip this loop except for domain 0
  if (listCounts != NULL) {
    for (dom = 0; dom < numDomains; dom++) {
      // see if domain 0 will send anything to domain dom
      if (listCounts[dom] > 0) {
        localMsgCnt[dom] = 1;
        localMsgCnt[numDomains] = 1;
      }
    }
  }

  DSMPI::AllReduce(localMsgCnt, globalMsgCnt, numDomains + 1, MPI_INT, MPI_SUM);

  // Do next send/receive of nodal data
  SendInitialNodeData(home, inData, globalMsgCnt, nodeLists, listCounts,
                      &nextAvailableTag);
  FreeInNodeArray(inData, param->nodeCount);
  FreeNodeLists(home, &nodeLists, &listCounts);

  /*
   *      This is a good place for a quick sanity check that the sum of
   *      nodes on all domains equals the total node count from the
   *      data file.
   */
  localNodeCount = 0;
  globalNodeCount = 0;

  for (i = 0; i < home->newNodeKeyPtr; i++) {
    if (home->nodeKeys[i] != (Node_t *)NULL) {
      localNodeCount++;
    }
  }

  DSMPI::Reduce(&localNodeCount, &globalNodeCount, 1, MPI_INT, MPI_SUM, 0);
  // Fatal("Stopped at reading file: ReadRestart.cpp"); //stop here

  if ((home->myDomain == 0) && (param->nodeCount != globalNodeCount)) {
    Fatal("ReadNodedataFile: Read %d nodes, expected %d!", globalNodeCount,
          param->nodeCount);
  }
  free(localMsgCnt);
  free(globalMsgCnt);
}

/* Function to read thermal stress
  --Yejun*/
void ReadStressDataFile(Home_t *home, InData_t *inData,
                        const char *stressFile) {
  int dom, iNbr;
  int maxTokenLen, tokenType, pIndex, okay;
  int valType, numVals;
  int numDomains;
  int numNbrs;
  int fileSeqNum = 0;
  int *globalMsgCnt, *localMsgCnt, **nodeLists, *listCounts;
  int localNodeCount, globalNodeCount;
  int binRead = 0;
  int nextAvailableTag = 0;
  real8 burgSumX, burgSumY, burgSumZ;
  void *valList;
  char inLine[500], token[256];
  char baseFileName[256], tmpFileName[256];
  FILE *fpSeg;
  Thermalstress_t *stress;
  Param_t *param;
  ParamList_t *stressParamList;

  param = home->param;
  numDomains = home->numDomains;
  fpSeg = (FILE *)NULL;

  stressParamList = home->stressParamList;
  globalMsgCnt = (int *)malloc((numDomains + 1) * sizeof(int));
  localMsgCnt = (int *)malloc((numDomains + 1) * sizeof(int));

  memset(inLine, 0, sizeof(inLine));
  maxTokenLen = sizeof(token);

  /*
   *      Only domain zero reads the initial stuff...
   */
  if (home->myDomain == 0) {
    if (stressFile == NULL) {
      Fatal("ReadStressDataFile(): ReadStressDataFile: No data file provided");
    }

    /*
     *          Try to open the nodal data file.  If the specified file
     *          can not be opened, try looking for a segmented data file
     *          (i.e. look for the first file segment <dataFile>.0). If
     *          that can't be opened either, then exit with an error.
     */

    if ((fpSeg = fopen(stressFile, "r")) == NULL) {
      Fatal(
          "ReadStressDataFile(): Error %d opening file %s to read stress data",
          errno, stressFile);
    } else {
      printf("Reading stress file: %s\n", stressFile);
    }

    /*
     *          Get the first token.  This should either be a known
     *          parameter identifier.
     */
    tokenType = GetNextToken(fpSeg, token, maxTokenLen);
    pIndex = LookupParam(stressParamList, token);

    /*
     *              Just go through the stress data file reading all
     *              the associated parameters.  Need to do special
     *              processing of the domain decomposition, and when
     *              when we hit the 'stressdata' identifier, just break
     *              out of this loop.
     */
    while ((tokenType != TOKEN_ERR) && (tokenType != TOKEN_NULL)) {
      if (pIndex >= 0) {
        /*
         *                      Token represents a known parameter identifier,
         * so read the associated value(s).
         */
        valType = stressParamList->varList[pIndex].valType;
        numVals = stressParamList->varList[pIndex].valCnt;
        valList = stressParamList->varList[pIndex].valList;
        okay = GetParamVals(fpSeg, valType, numVals, valList);
        if (!okay) {
          Fatal("ReadStressDataFile(): Parsing Error obtaining values for "
                "parameter %s\n",
                stressParamList->varList[pIndex].varName);
        }
      } else {
        /*
         *                      Token does not represent one of the simple
         *                      parameters.  If it's not one of the identifiers
         *                      that needs special handling, skip it.
         */
        if (strcmp(token, "stressData") == 0) {
          /*
           *                          When we hit the nodal data, we can break
           *                          out of the loop because we are assuming
           *                          all other data file parameters have been
           *                          processed.  If they have not, we have a
           *                          problem since processing of the nodal data
           *                          requires the other parameters.
           *                          Note: Remainder of the file should just
           *                          contain " = " followed by the nodal data,
           *                          so be sure to skip the next token before
           *                          processing the nodal data.
           */
          tokenType = GetNextToken(fpSeg, token, maxTokenLen);
          break;
        }
      }
      tokenType = GetNextToken(fpSeg, token, maxTokenLen);
      // printf("%d %s\n", tokenType,stressParamList->varList[pIndex].varName);
      if ((tokenType == TOKEN_NULL) || (tokenType == TOKEN_ERR)) {

        Fatal("ReadStressDataFile(): Parsing error on file %s\n", stressFile);
      }
      pIndex = LookupParam(stressParamList, token);
    }

    if (param->Time_evo < 1) {
      Fatal("Stress file generation format ERROR");
    }

  } /* if (home->myDomain == 0) */

  /*
   *      Domain zero now needs to pass various pieces of data to
   *      the remote domains before any processes begin to read
   *      and distribute the nodal data... start with the param structure.
   */
  DSMPI::Broadcast((char *)param, sizeof(Param_t), MPI_CHAR, 0);
  DSMPI::Barrier();
  // printf("it is %d %d %d %d %d %lf %lf %lf\n",home->myDomain,
  // param->Time_evo, param->GP_x, param->GP_y, param->GP_z,
  // param->stress_Dim[0], param->stress_Dim[1], param->stress_Dim[2]);

  unsigned int i, j, k, kk = 0;
  home->stress =
      (Thermalstress_t ****)calloc(param->GP_x, sizeof(Thermalstress_t ***));
  for (i = 0; i < param->GP_x; i++) {
    home->stress[i] =
        (Thermalstress_t ***)calloc(param->GP_y, sizeof(Thermalstress_t **));
    for (j = 0; j < param->GP_y; j++) {
      home->stress[i][j] =
          (Thermalstress_t **)calloc(param->GP_z, sizeof(Thermalstress_t *));
      for (k = 0; k < param->GP_z; k++) {
        home->stress[i][j][k] =
            (Thermalstress_t *)calloc(param->Time_evo, sizeof(Thermalstress_t));
      }
    }
  }

  int *home_stress_list_int;
  double *home_stress_list_double;
  home_stress_list_int = (int *)calloc(param->GP_x * param->GP_y * param->GP_z *
                                           param->Time_evo * 3,
                                       sizeof(int));
  home_stress_list_double = (double *)calloc(
      param->GP_x * param->GP_y * param->GP_z * param->Time_evo * 12,
      sizeof(double));

  if (home->myDomain == 0) {
    for (i = 0; i < param->GP_x; i++) {
      for (j = 0; j < param->GP_y; j++) {
        for (k = 0; k < param->GP_z; k++) {
          for (kk = 0; kk < param->Time_evo; kk++) {
            Getline(inLine, sizeof(inLine), fpSeg);
            if ((inLine[0] == 0) &&
                (i * param->GP_y * param->GP_z * param->Time_evo +
                     j * param->GP_z * param->Time_evo + k * param->Time_evo +
                     kk <
                 param->GP_x * param->GP_y * param->GP_z * param->Time_evo -
                     1)) {
              Fatal("ReadStressDataFile(): stress file ended before reading "
                    "all grid points\n");
            }

            stress = &home->stress[i][j][k][kk];
            sscanf(
                inLine,
                "%d %d %d, %lf %lf %lf, %lf %lf %lf %lf %lf %lf, %lf %lf, %lf",
                &stress->n_x, &stress->n_y, &stress->n_z, &stress->x,
                &stress->y, &stress->z, &stress->s_xx, &stress->s_yy,
                &stress->s_zz, &stress->s_yz, &stress->s_xz, &stress->s_xy,
                &stress->t_start, &stress->t_end, &stress->Temp);
            // printf("%d %d %d %d %d %d
            // %d\n",stress->n_x,stress->n_y,stress->n_z,i,j,k,kk);
            if (stress->n_x != i + 1 || stress->n_y != j + 1 ||
                stress->n_z != k + 1) {
              // printf("%d %d %d, %d %d %d, %d %lf %lf %lf %lf \n",
              // stress->n_x, stress->n_y, stress->n_z, i, j, k, kk, stress->x,
              // stress->y, stress->z, stress->t_end);
              Fatal("Stress file mesh generation format ERROR");
            }

            // make sure it is the primary image
            GetPrimaryImage(home->param, &stress->x, &stress->y, &stress->z);

            home_stress_list_int[i * param->GP_y * param->GP_z *
                                     param->Time_evo +
                                 j * param->GP_z * param->Time_evo +
                                 k * param->Time_evo + kk] = stress->n_x;
            home_stress_list_int
                [i * param->GP_y * param->GP_z * param->Time_evo +
                 j * param->GP_z * param->Time_evo + k * param->Time_evo + kk +
                 param->GP_x * param->GP_y * param->GP_z * param->Time_evo] =
                    stress->n_y;
            home_stress_list_int[i * param->GP_y * param->GP_z *
                                     param->Time_evo +
                                 j * param->GP_z * param->Time_evo +
                                 k * param->Time_evo + kk +
                                 param->GP_x * param->GP_y * param->GP_z *
                                     param->Time_evo * 2] = stress->n_z;

            home_stress_list_double[i * param->GP_y * param->GP_z *
                                        param->Time_evo +
                                    j * param->GP_z * param->Time_evo +
                                    k * param->Time_evo + kk] = stress->x;
            home_stress_list_double
                [i * param->GP_y * param->GP_z * param->Time_evo +
                 j * param->GP_z * param->Time_evo + k * param->Time_evo + kk +
                 param->GP_x * param->GP_y * param->GP_z * param->Time_evo] =
                    stress->y;
            home_stress_list_double[i * param->GP_y * param->GP_z *
                                        param->Time_evo +
                                    j * param->GP_z * param->Time_evo +
                                    k * param->Time_evo + kk +
                                    param->GP_x * param->GP_y * param->GP_z *
                                        param->Time_evo * 2] = stress->z;
            home_stress_list_double[i * param->GP_y * param->GP_z *
                                        param->Time_evo +
                                    j * param->GP_z * param->Time_evo +
                                    k * param->Time_evo + kk +
                                    param->GP_x * param->GP_y * param->GP_z *
                                        param->Time_evo * 3] = stress->s_xx;
            home_stress_list_double[i * param->GP_y * param->GP_z *
                                        param->Time_evo +
                                    j * param->GP_z * param->Time_evo +
                                    k * param->Time_evo + kk +
                                    param->GP_x * param->GP_y * param->GP_z *
                                        param->Time_evo * 4] = stress->s_yy;
            home_stress_list_double[i * param->GP_y * param->GP_z *
                                        param->Time_evo +
                                    j * param->GP_z * param->Time_evo +
                                    k * param->Time_evo + kk +
                                    param->GP_x * param->GP_y * param->GP_z *
                                        param->Time_evo * 5] = stress->s_zz;
            home_stress_list_double[i * param->GP_y * param->GP_z *
                                        param->Time_evo +
                                    j * param->GP_z * param->Time_evo +
                                    k * param->Time_evo + kk +
                                    param->GP_x * param->GP_y * param->GP_z *
                                        param->Time_evo * 6] = stress->s_yz;
            home_stress_list_double[i * param->GP_y * param->GP_z *
                                        param->Time_evo +
                                    j * param->GP_z * param->Time_evo +
                                    k * param->Time_evo + kk +
                                    param->GP_x * param->GP_y * param->GP_z *
                                        param->Time_evo * 7] = stress->s_xz;
            home_stress_list_double[i * param->GP_y * param->GP_z *
                                        param->Time_evo +
                                    j * param->GP_z * param->Time_evo +
                                    k * param->Time_evo + kk +
                                    param->GP_x * param->GP_y * param->GP_z *
                                        param->Time_evo * 8] = stress->s_xy;
            home_stress_list_double[i * param->GP_y * param->GP_z *
                                        param->Time_evo +
                                    j * param->GP_z * param->Time_evo +
                                    k * param->Time_evo + kk +
                                    param->GP_x * param->GP_y * param->GP_z *
                                        param->Time_evo * 9] = stress->t_start;
            home_stress_list_double[i * param->GP_y * param->GP_z *
                                        param->Time_evo +
                                    j * param->GP_z * param->Time_evo +
                                    k * param->Time_evo + kk +
                                    param->GP_x * param->GP_y * param->GP_z *
                                        param->Time_evo * 10] = stress->t_end;
            home_stress_list_double[i * param->GP_y * param->GP_z *
                                        param->Time_evo +
                                    j * param->GP_z * param->Time_evo +
                                    k * param->Time_evo + kk +
                                    param->GP_x * param->GP_y * param->GP_z *
                                        param->Time_evo * 11] = stress->Temp;
            // printf("it is %d %d %d %lf %lf %lf %lf\n",stress->n_x,
            // stress->n_y, stress->n_z, stress->x, stress->y, stress->z,
            // stress->s_xy);
          }
        }
      }
    }
  }
  DSMPI::Barrier();
  DSMPI::Broadcast(home_stress_list_int, param->GP_x * param->GP_y *
                                             param->GP_z * param->Time_evo * 3,
                   MPI_INT, 0);
  DSMPI::Broadcast(home_stress_list_double,
                   param->GP_x * param->GP_y * param->GP_z * param->Time_evo *
                       12,
                   MPI_DOUBLE, 0);
  DSMPI::Barrier();
  if (home->myDomain != 0) {
    for (i = 0; i < param->GP_x; i++) {
      for (j = 0; j < param->GP_y; j++) {
        for (k = 0; k < param->GP_z; k++) {
          for (kk = 0; kk < param->Time_evo; kk++) {

            home->stress[i][j][k][kk].n_x =
                home_stress_list_int[i * param->GP_y * param->GP_z *
                                         param->Time_evo +
                                     j * param->GP_z * param->Time_evo +
                                     k * param->Time_evo + kk];
            home->stress[i][j][k][kk].n_y = home_stress_list_int
                [i * param->GP_y * param->GP_z * param->Time_evo +
                 j * param->GP_z * param->Time_evo + k * param->Time_evo + kk +
                 param->GP_x * param->GP_y * param->GP_z * param->Time_evo];
            home->stress[i][j][k][kk].n_z = home_stress_list_int
                [i * param->GP_y * param->GP_z * param->Time_evo +
                 j * param->GP_z * param->Time_evo + k * param->Time_evo + kk +
                 param->GP_x * param->GP_y * param->GP_z * param->Time_evo * 2];

            home->stress[i][j][k][kk].x =
                home_stress_list_double[i * param->GP_y * param->GP_z *
                                            param->Time_evo +
                                        j * param->GP_z * param->Time_evo +
                                        k * param->Time_evo + kk];
            home->stress[i][j][k][kk].y = home_stress_list_double
                [i * param->GP_y * param->GP_z * param->Time_evo +
                 j * param->GP_z * param->Time_evo + k * param->Time_evo + kk +
                 param->GP_x * param->GP_y * param->GP_z * param->Time_evo];
            home->stress[i][j][k][kk].z = home_stress_list_double
                [i * param->GP_y * param->GP_z * param->Time_evo +
                 j * param->GP_z * param->Time_evo + k * param->Time_evo + kk +
                 param->GP_x * param->GP_y * param->GP_z * param->Time_evo * 2];

            home->stress[i][j][k][kk].s_xx = home_stress_list_double
                [i * param->GP_y * param->GP_z * param->Time_evo +
                 j * param->GP_z * param->Time_evo + k * param->Time_evo + kk +
                 param->GP_x * param->GP_y * param->GP_z * param->Time_evo * 3];
            home->stress[i][j][k][kk].s_yy = home_stress_list_double
                [i * param->GP_y * param->GP_z * param->Time_evo +
                 j * param->GP_z * param->Time_evo + k * param->Time_evo + kk +
                 param->GP_x * param->GP_y * param->GP_z * param->Time_evo * 4];
            home->stress[i][j][k][kk].s_zz = home_stress_list_double
                [i * param->GP_y * param->GP_z * param->Time_evo +
                 j * param->GP_z * param->Time_evo + k * param->Time_evo + kk +
                 param->GP_x * param->GP_y * param->GP_z * param->Time_evo * 5];
            home->stress[i][j][k][kk].s_yz = home_stress_list_double
                [i * param->GP_y * param->GP_z * param->Time_evo +
                 j * param->GP_z * param->Time_evo + k * param->Time_evo + kk +
                 param->GP_x * param->GP_y * param->GP_z * param->Time_evo * 6];
            home->stress[i][j][k][kk].s_xz = home_stress_list_double
                [i * param->GP_y * param->GP_z * param->Time_evo +
                 j * param->GP_z * param->Time_evo + k * param->Time_evo + kk +
                 param->GP_x * param->GP_y * param->GP_z * param->Time_evo * 7];
            home->stress[i][j][k][kk].s_xy = home_stress_list_double
                [i * param->GP_y * param->GP_z * param->Time_evo +
                 j * param->GP_z * param->Time_evo + k * param->Time_evo + kk +
                 param->GP_x * param->GP_y * param->GP_z * param->Time_evo * 8];
            home->stress[i][j][k][kk].t_start = home_stress_list_double
                [i * param->GP_y * param->GP_z * param->Time_evo +
                 j * param->GP_z * param->Time_evo + k * param->Time_evo + kk +
                 param->GP_x * param->GP_y * param->GP_z * param->Time_evo * 9];
            home->stress[i][j][k][kk].t_end =
                home_stress_list_double[i * param->GP_y * param->GP_z *
                                            param->Time_evo +
                                        j * param->GP_z * param->Time_evo +
                                        k * param->Time_evo + kk +
                                        param->GP_x * param->GP_y *
                                            param->GP_z * param->Time_evo * 10];
            home->stress[i][j][k][kk].Temp =
                home_stress_list_double[i * param->GP_y * param->GP_z *
                                            param->Time_evo +
                                        j * param->GP_z * param->Time_evo +
                                        k * param->Time_evo + kk +
                                        param->GP_x * param->GP_y *
                                            param->GP_z * param->Time_evo * 11];
          }
        }
      }
    }
  }

  free(home_stress_list_int);
  free(home_stress_list_double);
  home_stress_list_int = NULL;
  home_stress_list_double = NULL;
}
