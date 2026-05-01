/****************************************************************************
 *
 *      Module:      Decomp.c
 *      Description: Contains a set of generic functions used to create,
 *                   access, and manipulate the domain decomposition
 *                   regardless of the type of decomposition in use.
 *                   These generic functions will do any necessary
 *                   setup and invoke the appropriate functions to do
 *                   the real work based on the type of domain decomposition
 *                   currently in use.
 *
 *      Includes public functions
 *          BroadcastDecomp()
 *          FindCoordDomain()
 *          FreeDecomp()
 *          GetAllDecompBounds()
 *          GetCellDomainList()
 *          GetLocalDomainBounds()
 *          ReadDecompBounds()
 *          Rebalance()
 *          UniformDecomp()
 *          WriteDecompBounds()
 *
 *      Includes private functions
 *          DLBStats()
 *          GetLoadData()
 *
 ***************************************************************************/
#include "Home.h"
#include "Restart.h"
#include "Decomp.h"
#include "RSDecomp.h"
#include "DSMPI.h"

/*
 *      Define some flags that toggle compile-time generation of
 *      some code for debugging the load-balancing algorithms.
 *      For each definition, 0 == OFF, 1 == ON.
 *
 *      DLB_DEBUG           Highest level flag enabling general debug code.
 *      DLB_PRINT_STATS     Enables message printed to stdout each DLB cycle
 *                          with load-balance stats such as min and max
 *                          per-process load, average load, and so on.  This
 *                          has no effect when DLB_DEBUG == 0.
 *      DLB_DUMP_DATAFILES  Enables creation of a file
 *                          (DLB_stats.<sequnce_num>) every DLB cycle with
 *                          the per-process load data.  This has no effect
 *                          when DLB_DEBUG == 0.
 */
#define DLB_DEBUG 0
#define DLB_PRINT_STATS 1
#define DLB_DUMP_DATAFILES 0

void NodeBasedBalance(Home_t *poHome);

/*-----------------------------------------------------------------------
 *
 *      Function:     DLBStats
 *      Description:  This function is only provided for debugging
 *                    purposes.  It handles writing (to stdout and/or
 *                    disk) some data related to the current load
 *                    on all processors.  Only domain zero will
 *                    do any output here...
 *
 *      Arguments:
 *          loadData  Array containing current load data (wallclock
 *                    time or force calc counts) for each domain.
 *                    Array is assumed to have same number of elements
 *                    as the job's domain count.
 *
 *----------------------------------------------------------------------*/
static void DLBStats(Home_t *home, real8 *loadData) {
  int domain;
  real8 minVal, maxVal, totVal, avgVal, imbalance;
  char statFileName[64];
  FILE *statFP;
  static int seqNum = 1;

  if (home->myDomain != 0) {
    return;
  }

#if DLB_DUMP_DATAFILES
  /*
   *      If enabled, dump the current load data for all processors
   *      to a disk file.
   */
  sprintf(statFileName, "DLB_stats.%d", seqNum);
  statFP = fopen(statFileName, "w");
  if (statFP == (FILE *)NULL) {
    printf("  *** Error %d opening file %s\n", errno, statFileName);
    return;
  }

  for (domain = 0; domain < home->numDomains; domain++) {
    fprintf(statFP, "%d  %lf\n", domain, loadData[domain]);
  }

  fclose(statFP);
#endif

#if DLB_PRINT_STATS
  /*
   *      If enabled, print to stdout a single line with some stats
   *      about the current load balance.
   */
  minVal = loadData[0];
  maxVal = 0.0;
  totVal = 0.0;
  avgVal = 0.0;

  for (domain = 0; domain < home->numDomains; domain++) {
    if (loadData[domain] < minVal)
      minVal = loadData[domain];
    if (loadData[domain] > maxVal)
      maxVal = loadData[domain];
    totVal += loadData[domain];
  }

  avgVal = totVal / home->numDomains;

  if (avgVal == 0.0) {
    imbalance = 0.0;
  } else {
    imbalance = (maxVal - avgVal) / avgVal;
  }

  printf(" +++ DLB %d: min=%lf max=%lf avg=%lf imbalance=%lf\n", seqNum, minVal,
         maxVal, avgVal, imbalance);
#endif

  seqNum++;

  return;
}

/*------------------------------------------------------------------------
 *
 *      Function:       BroadcastDecomp
 *      Description:    Generic function to control broadcast of the
 *                      domain decomposition from task zero to all other
 *                      tasks.  This function will determine the type of
 *                      domain decomposition is use and invoke the
 *                      broadcast function appropriate to the decomposition.
 *
 *      Arguments:
 *          decomp   Pointer to decomposition to be broadcast from domain
 *                   zero to all other domains. Pointer is NULL or
 *                   uninitialized on all tasks but zero.
 *
 *-----------------------------------------------------------------------*/
void BroadcastDecomp(Home_t *home, void *decomp) {
  BroadcastRSDecomp(home, (RSDecomp_t *)decomp);
  /*
   *      Have each process copy its own local boundaries into the
   *      <home> structure... saves multiple lookups later on.
   */
  GetLocalDomainBounds(home, home->decomp);
}

/*---------------------------------------------------------------------------
 *
 *      Function:       UniformRBDecomp
 *      Description:    Generic function to generate a uniform domain
 *                      decomposition in which all domains will encompass
 *                      volumes of equivalent size and shape.  This function
 *                      will determine the type of domain decomposition in
 *                      use and invoke the appropriate decomposition function.
 *
 *      Arguments:
 *          decomp  Location in which to return to the caller a pointer
 *                  to a new uniform domain decomposition
 *
 *-------------------------------------------------------------------------*/
void UniformDecomp(Home_t *poHome, void **decomp) {
  UniformRSDecomp(poHome->param, (RSDecomp_t **)decomp);
}

RSDecomp_t *PointBasedDecomp(Home_t *poHome, InData_t *poInputData) {
  unsigned int iXDomainsCount = poInputData->param->nXdoms;
  unsigned int iYDomainsCount = poInputData->param->nYdoms;
  unsigned int iZDomainsCount = poInputData->param->nZdoms;
  list<AxisAlignedBoundingBox *> lpoDomains;
  double dXMin = -0.5 * poInputData->param->Dimensions[0];
  double dYMin = -0.5 * poInputData->param->Dimensions[1];
  double dZMin = -0.5 * poInputData->param->Dimensions[2];
  double dXMax = 0.5 * poInputData->param->Dimensions[0];
  double dYMax = 0.5 * poInputData->param->Dimensions[1];
  double dZMax = 0.5 * poInputData->param->Dimensions[2];
  AxisAlignedBoundingBox oBox;
  oBox.SetXMin(dXMin);
  oBox.SetYMin(dYMin);
  oBox.SetZMin(dZMin);
  oBox.SetXMax(dXMax);
  oBox.SetYMax(dYMax);
  oBox.SetZMax(dZMax);
  list<Point> loPoints;
  Node_t *poNode = NULL;
  unsigned int i = 0;
  for (i = 0; i < poInputData->param->nodeCount; i++) {
    poNode = &poInputData->node[i];
    loPoints.push_back(GetInLocalCoordinates(
        poHome->param, Point(poNode->x, poNode->y, poNode->z)));
  }
  lpoDomains = oBox.PointBasedPartition(&loPoints, iXDomainsCount,
                                        iYDomainsCount, iZDomainsCount);

  loPoints.clear();
  unsigned int j = 0;
  unsigned int k = 0;
  double *newXBounds = (double *)malloc((iXDomainsCount + 1) * sizeof(double));
  double **newYBounds = (double **)malloc(iXDomainsCount * sizeof(double *));
  double ***newZBounds = (double ***)malloc(iXDomainsCount * sizeof(double **));
  list<AxisAlignedBoundingBox *>::iterator liDomains = lpoDomains.begin();

  for (i = 0; i < iXDomainsCount; i++) {
    newXBounds[i] = (*liDomains)->GetXMin();
    newYBounds[i] = (double *)malloc((iYDomainsCount + 1) * sizeof(double));
    newZBounds[i] = (double **)malloc(iYDomainsCount * sizeof(double *));
    for (j = 0; j < iYDomainsCount; j++) {
      newYBounds[i][j] = (*liDomains)->GetYMin();
      newZBounds[i][j] =
          (double *)malloc((iZDomainsCount + 1) * sizeof(double));
      for (k = 0; k < iZDomainsCount; k++) {
        newZBounds[i][j][k] = (*liDomains)->GetZMin();
        delete (*liDomains);
        liDomains++;
      }
      newZBounds[i][j][k] = dZMax;
    }
    newYBounds[i][j] = dYMax;
  }
  newXBounds[i] = dXMax;
  lpoDomains.clear();
  RSDecomp_t *decomp = (RSDecomp_t *)malloc(sizeof(RSDecomp_t));
  decomp->domBoundX = newXBounds;
  decomp->domBoundY = newYBounds;
  decomp->domBoundZ = newZBounds;
  return decomp;
}
/*-------------------------------------------------------------------------
 *
 *      Function:    WriteDecompBounds
 *      Description: Generic function to write the domain decomposition
 *                   into a restart file.  This function will determine
 *                   the type of domain decomposition in use and invoke
 *                   the write function appropriate to the decomposition
 *                   type.
 *
 *      Arguments:
 *          fp          open file descriptor for the restart file being
 *                      written
 *          versionNum  version number of the nodal data file being written.
 *
 *------------------------------------------------------------------------*/
void WriteDecompBounds(Home_t *home, FILE *fp) {
  WriteRSDecompBounds(home, fp, (RSDecomp_t *)home->decomp);
}

/*-------------------------------------------------------------------------
 *
 *      Function:    GetAllDecompBounds
 *      Description: Generic function to create a single dimension array
 *                   containing all the domain boundaries.  This function
 *                   will determine the type of domain decomposition in
 *                   use, calculate the necessary array size and invoke
 *                   the 'get' function appropriate to the decomposition
 *                   type.
 *
 *      Arguments:
 *          bounds      location in whcih to return to the caller the
 *                      pointer to the array of domain boundaries.
 *          numBounds   location in which to return to the caller the
 *                      number of elements in the <*bounds> array.
 *
 *------------------------------------------------------------------------*/
void GetAllDecompBounds(Home_t *home, real8 **bounds, int *numBounds) {
  int xDoms, yDoms, zDoms;
  xDoms = home->param->nXdoms;
  yDoms = home->param->nYdoms;
  zDoms = home->param->nZdoms;
  *numBounds = xDoms * (yDoms * (zDoms + 1) + (yDoms + 1)) + (xDoms + 1);
  *bounds = (real8 *)malloc(*numBounds * sizeof(real8));
  GetAllRSDecompBounds(home, (RSDecomp_t *)home->decomp, *bounds);
}

/*-------------------------------------------------------------------------
 *
 *      Function:    FreeDecomp
 *      Description: Generic function to release memory associated with
 *                   a domain decomposition.  This function will
 *                   invoke the memory freeing function appropriate
 *                   to the decomposition type.
 *
 *      Arguments:
 *          decomp   pointer to the domain decomposition to be freed
 *
 *------------------------------------------------------------------------*/
void FreeDecomp(Home_t *home, void *decomp) {
  FreeRSDecomp(home, (RSDecomp_t *)decomp);
  free(decomp);
}

/*-------------------------------------------------------------------------
 *
 *      Function:    GetLocalDomainBounds
 *      Description: Generic function to search the domain decomposition
 *                   for the current domain's boundaries and store them
 *                   in the <home> structure (to save multiple lookups
 *                   later).
 *
 *      Arguments:
 *          decomp   pointer to the current domain decomposition
 *
 *------------------------------------------------------------------------*/
void GetLocalDomainBounds(Home_t *home, void *decomp) {
  GetRSDecompLocalDomainBounds(home, (RSDecomp_t *)decomp);
}

/*-------------------------------------------------------------------------
 *
 *      Function:    GetCellDomainList
 *      Description: Generic function to search the domain decomposition
 *                   for the list of domains intersecting the specified
 *                   cell.  This function will invoke the search function
 *                   appropriate to the type of domain decomposition
 *                   being used.
 *
 *      Arguments:
 *          cellID    ID of cell as returned by EncodeCellIdx().
 *          domCount  Location in which to return to caller the number of
 *                    domains intersecting the specified cell.
 *          domList   Location in which to return to caller the array
 *                    containing the IDs of all domains intersecting the
 *                    specified cell.
 *
 *------------------------------------------------------------------------*/
void GetCellDomainList(Home_t *home, int cellID, int *domCount, int **domList) {
  GetRSDecompCellDomainList(home, cellID, domCount, domList);
}

/*-------------------------------------------------------------------------
 *
 *      Function:    GetLoadData
 *      Description: Obtains the per-process load information required
 *                   to do dynamic load-balancing across the processors.
 *
 *      Arguments:
 *          criteria  Identifies the type of load data to obtain.
 *                    DLB_USE_FORCECALC_COUNT = number of force calculations
 *                                              done this cycle
 *                    DLB_USE_WALLCLK_TIME    = wallclock time spent doing
 *                                              force calculations this
 *                                              cycle.
 *          loadData  Location in which to return to the caller a pointer
 *                    to the load data array.  The array consists of 1
 *                    value per process.
 *
 *------------------------------------------------------------------------*/
static void GetLoadData(Home_t *home, real8 **loadData) {
  real8 localLoad = (real8)home->cycleForceCalcCount;
  if (localLoad < 1.0) {
    localLoad = 1.0;
  }
  real8 *globalLoadData;
  /*
   *      Caller is responsible for freeing this buffer!
   */
  globalLoadData = (real8 *)malloc(home->numDomains * sizeof(real8));

#ifdef PARALLEL
  MPI_Allgather(&localLoad, 1, MPI_DOUBLE, globalLoadData, 1, MPI_DOUBLE,
                MPI_COMM_WORLD);
#else
  *globalLoadData = localLoad;
#endif
  *loadData = globalLoadData;
}

/*-------------------------------------------------------------------------
 *
 *      Function:    Rebalance
 *      Description: Generic function to dynamically rebalance the
 *                   work load by adjusting the domain decomposition.
 *                   This function will obtain updated per-process load
 *                   data which it will then pass on to the decomposition
 *                   (and related rebalancing) functions appropriate to
 *                   the type of domain decomposition being used.
 *
 *      Arguments:
 *          criteria  Identifies the type of per process load data to obtain.
 *                    DLB_USE_FORCECALC_COUNT = number of force calculations
 *                                              done this cycle
 *                    DLB_USE_WALLCLK_TIME    = wallclock time spent doing
 *                                              force calculations this
 *                                              cycle.
 *------------------------------------------------------------------------*/
void Rebalance(Home_t *home) {
  // NodeBasedBalance(home);
  int didDLB = 0;
  real8 *loadData;
  Param_t *param;

  param = home->param;

  /*
   *      There's a number of conditions that determine whether to do
   *      load balancing this cycle... check 'em all.
   */
  if ((param->DLBfreq > 0) && (home->cycle % param->DLBfreq <= 3)) {
    /*
     *          Decomposition requires the raw per-process load
     *          data, so get the data (and dump statistics about
     *          the current load imbalance if necessary).
     */
    GetLoadData(home, &loadData);

    /*
     *              For recursive sectioning domain decompositions, load
     *              balancing occurs over a 3-cycle period; balance along
     *              X dimension the first cycle, along the Y dimension the
     *              second cycle and along Z the third cycle.
     */
    switch (home->cycle % param->DLBfreq) {
    case 0:
      DLBalanceX(home, loadData);
      break;
    case 1:
      DLBalanceY(home, loadData);
      break;
    case 2:
      DLBalanceZ(home, loadData);
      break;
    }

    GetLocalDomainBounds(home, home->decomp);
    didDLB = 1;

    /*
     *          If a rebalance was actually performed, then free
     *          the old cell and neighbor domain structures (since
     *          cell membership, neighbors, and neighboring domains are
     *          dependent on the domain boundaries) and reinitialize
     *          the structures.
     */
    free(loadData);

    if (didDLB) {
      DLBfreeOld(home);
      InitCellNatives(home);
      InitCellNeighbors(home);
      InitCellDomains(home);
      InitRemoteDomains(home);
    }
    // finally, update APB cell ownership if applicable
    home->poPrecipitateServer->UpdateAPBCellOwnership(home);
  }
}

/*-------------------------------------------------------------------------
 *
 *      Function:    FindCoordDomain
 *      Description: Generic function to search the domain decomposition
 *                   to locate the domain 'owning' the specified
 *                   coordinates.  This function will shift the specified
 *                   coordinates into the primary image space (if
 *                   necessary) then invoke the search function appropriate
 *                   to the type of domain decomposition being used.
 *
 *      Arguments:
 *          updateCoords Flag indicating if coordinates outside the primary
 *                       image should be converted to the corresponding
 *                       coordinates.  Any non-zero value permits the
 *                       coordinate update.  Note: regardless of the
 *                       value of this flag, no coordinate update will
 *                       be done if periodic boundaries are not enabled.
 *          x, y, z      Pointers to coordinates.  Contents may be updated
 *                       depending on the value of <updateCoords>.
 *
 *------------------------------------------------------------------------*/
int FindCoordDomain(Home_t *home, const double &x, const double &y,
                    const double &z) {
  int domID;
  real8 newX, newY, newZ;
  real8 xMin, yMin, zMin;
  real8 xMax, yMax, zMax;
  Param_t *param;

  param = home->param;
  domID = home->myDomain;

  newX = x;
  newY = y;
  newZ = z;

  GetInLocalCoordinates(home->param, newX, newY, newZ);

  /*
   *      Find the min/max coordinates of the current domain.  If the
   *      specified position is contained within the current domain
   *      we don't need to go any further.
   */
  xMin = home->domXmin;
  yMin = home->domYmin;
  zMin = home->domZmin;

  xMax = home->domXmax;
  yMax = home->domYmax;
  zMax = home->domZmax;

  if (((newX >= xMin) && (newX < xMax)) && ((newY >= yMin) && (newY < yMax)) &&
      ((newZ >= zMin) && (newZ < zMax))) {
    return (domID);
  }

  domID = FindRSDecompCoordDomain(home, (RSDecomp_t *)home->decomp, newX, newY,
                                  newZ);
  return domID;
}
int FindCoordDomain(Home_t *home, Point *poPoint) {
  return FindCoordDomain(home, poPoint->GetX(), poPoint->GetY(),
                         poPoint->GetZ());
}

void NodeBasedBalance(Home_t *poHome) {
  DSMPI::Barrier();
  if (poHome->myDomain == 0) {
    char cLatestRestartFileName[512];
    sprintf(cLatestRestartFileName, "%s/rs%04d%s", DIR_RESTART,
            poHome->param->savecounter, NODEDATA_FILE_SUFFIX);
    FILE *fpLatestRestart = fopen(cLatestRestartFileName, "r");
    string sRead;
    unsigned int iNodeCount = 0;
    // skip the header section
    while (!feof(fpLatestRestart)) {
      sRead = GetRealString(1024, fpLatestRestart);
      if (sRead.compare(0, 9, "nodeCount", 9) == 0) {
        sscanf(sRead.c_str(), "nodeCount =%d\n", &iNodeCount);
      }
      if (sRead.compare(0, 9, "nodalData", 9) == 0) {
        break;
      }
    }
    sRead = GetRealString(1024, fpLatestRestart);
    vector<double> vdXCoordinates(iNodeCount);
    vector<double> vdYCoordinates(iNodeCount);
    vector<double> vdZCoordinates(iNodeCount);
    unsigned int i = 0;
    unsigned int j = 0;
    unsigned int iNeighboursCount = 0;
    for (i = 0; i < iNodeCount; i++) {
      sRead = GetRealString(1024, fpLatestRestart);
      sscanf(sRead.c_str(), "%*d,%*d\t%lf\t%lf\t%lf\t%d\t%*d\n",
             &vdXCoordinates[i], &vdXCoordinates[i], &vdXCoordinates[i],
             &iNeighboursCount);
      // skip surface normal
      sRead = GetRealString(1024, fpLatestRestart);
      for (j = 0; j < iNeighboursCount; j++) {
        // skip arm data
        sRead = GetRealString(1024, fpLatestRestart);
        sRead = GetRealString(1024, fpLatestRestart);
      }
    }
    SupportSystem::QuickSort<double>(vdXCoordinates);
    fclose(fpLatestRestart);
  }
  DSMPI::Barrier();
}
