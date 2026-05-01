/****************************************************************************
 *
 *      Author:       Gregg Hommes
 *
 *      Module:       Initialize.c
 *
 *      Description:  Contains the driver routine for the initialization
 *                    of the application. Handles some of the general
 *                    initializations directly and calls all the more
 *                    specific initialization routines.
 *
 *      Last Modified: 01/09/08: Gregg Hommes - Added VerifyBurgersVectors()
 *                               sanity check.
 *
 ****************************************************************************/

#include <memory.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/utsname.h>
#include <fcntl.h>
#include <ctype.h>
#include <time.h>
#include <pwd.h>
#include "Home.h"
#include "Init.h"
#include "InData.h"
#include "Mobility.h"
#include "Decomp.h"
#include "Parse.h"
#include "Restart.h"
#include "Comm.h"
#include "DSMPI.h"

/*---------------------------------------------------------------------------
 *
 *      Function:     InitRecycleNodeHeap
 *      Description:  If the user requested preservation (if possible)
 *                    of the node tags found in the restart file, the
 *                    <home->nodeKeys> array may be sparsely populated right
 *                    from the start.  In this case, we have to
 *                    create an initial heap of recycled nodes containing
 *                    the indices of all unused <home->nodeKeys> entries
 *                    less than <home->newNodeKeyPtr>
 *
 *-------------------------------------------------------------------------*/
static void InitRecycleNodeHeap(Home_t *home) {
  int i;

  for (i = 0; i < home->newNodeKeyPtr; i++) {
    if (home->nodeKeys[i] == (Node_t *)NULL) {
      RecycleNodeTag(home, i);
    }
  }

  return;
}

/*---------------------------------------------------------------------------
 *
 *      Function:     OpenDir
 *      Description:  This function will create (if they does not yet exist)
 *                    the primary output directory for the run, plus all
 *                    necessary subdirectories for specific output types.
 *
 *-------------------------------------------------------------------------*/
int OpenDir(Home_t *home) {
  char *dirname = home->param->dirname;
  char subdir[256];

  /*
   *      Only domain zero creates the primary output directory; don't
   *      want thousands of processors beating on the file system.
   */
  if (home->myDomain == 0) {
    if (mkdir(dirname, S_IRWXU) != 0) {
      if (errno == EEXIST) {
        printf("Warning: %s already exists\n", dirname);
      } else {
        Fatal("Open error %d on directory %s", errno, dirname);
      }
    }
  }

  /*
   *      All processes wait for task zero to create the primary output
   *      directory then cd into that directory.
   */

  DSMPI::Barrier();

  if (chdir(dirname) != 0) {
    Fatal("Task %d: Unable to cd into directory %s", home->myDomain, dirname);
  }

  DSMPI::Barrier();

  if (home->myDomain == 0)
    printf("chdir successful on all tasks.\n");

  /*
   *      Create all subdirectories needed for specific types of output.
   *      Again, only domain zero need do these creates.
   *
   *      Note: The current working directory for all tasks is now the
   *      user specified output directory, so when we create the
   *      subdirectories for various output types we just create them
   *      local to the current working directory.
   */
  if (home->myDomain == 0) {
    snprintf(subdir, sizeof(subdir), "./%s", DIR_PROPERTIES);
    (void)mkdir(subdir, S_IRWXU);

    snprintf(subdir, sizeof(subdir), "./%s", DIR_RESTART);
    (void)mkdir(subdir, S_IRWXU);

    if (home->param->tecplot) {
      snprintf(subdir, sizeof(subdir), "./%s", DIR_TECPLOT);
      (void)mkdir(subdir, S_IRWXU);
    }
  }

  return (0);
}

/*---------------------------------------------------------------------------
 *
 *      Function:     SetRemainingDefaults
 *      Description:  The default values of certain global parameters
 *                    are special in that they depend on values of
 *                    other global parameters.  If the user did not
 *                    specify values for these special parameters,
 *                    this function will calculate the necessary
 *                    defaults (as well as do some additional sanity
 *                    checks on some of the values).
 *
 *-------------------------------------------------------------------------*/
void SetRemainingDefaults(Home_t *home) {
  real8 tmp, eps;
  real8 xCellSize, yCellSize, zCellSize, minCellSize;
  Param_t *param;

  param = home->param;

  param->delSegLength = 0.0;

  xCellSize = param->Lx / param->nXcells;
  yCellSize = param->Ly / param->nYcells;
  zCellSize = param->Lz / param->nZcells;

  minCellSize = MIN(xCellSize, yCellSize);
  minCellSize = MIN(minCellSize, zCellSize);

  eps = 1.0e-02;

  /*
   *      The core radius and maximum segment length are required
   *      inputs.  If the user did not provide both values, abort
   *      now.
   */
  if (home->myDomain == 0) {
    if (param->rc < 0.0) {
      Fatal("The <rc> parameter is required but was not \n"
            "    provided in the control file");
    }

    if (param->maxSeg < 0.0) {
      Fatal("The <maxSeg> parameter is required but was not \n"
            "    provided in the control file");
    }
  }

  /*
   *      If not provided, set position error tolerance based on <rc>
   */
  if (param->rTol <= 0.0) {
    param->rTol = 0.25 * param->rc;
  }

  /*
   *      The deltaTT is set in the timestep integrator, but some
   *      mobility functions now use the deltaTT value, so it must
   *      be initialized before ParadisStep() is called since there
   *      is an initial mobility calculation done *before* timestep
   *      integration the first time into the function. Use this initial value
   */
  param->deltaTT = 1.0e-11;

  /*
   *      Set annihilation distance based on <rc>
   */
  param->rann = 2.0 * param->rTol;

  /*
   *      Minimum area criteria for remesh is dependent on maximum
   *      and minumum segment lengths and position error tolerance.
   */
  param->remeshAreaMin = 2.0 * param->rTol * param->maxSeg;

  if (param->minSeg > 0.0) {
    param->remeshAreaMin = MIN(param->remeshAreaMin,
                               (param->minSeg * param->minSeg * sqrt(3.0) / 4));
  }

  /*
   *      Maximum area criteria for remesh is dependent on minimum area,
   *      and maximum segment length.
   */
  param->remeshAreaMax =
      0.5 * ((4.0 * param->remeshAreaMin) +
             (0.25 * sqrt(3.0)) * (param->maxSeg * param->maxSeg));

  /*
   *      If the user did not provide a minSeg length, calculate one
   *      based on the remesh minimum area criteria.
   */
  if (param->minSeg <= 0.0) {
    param->minSeg = sqrt(param->remeshAreaMin * (4.0 / sqrt(3.0)));
  }

  /*
   *      If the user did not provide an Ecore value, set the default
   *      based on the shear modulus and rc values
   */
  if (param->Ecore < 0.0) {
    param->Ecore = (param->shearModulus / (4 * M_PI)) * log(param->rc / 0.1);
  }

  /*
   *      Now do some additional sanity checks.
   */
  if (home->myDomain == 0) {

    /*
     *          First check for some fatal errors...
     */
    if (param->maxSeg <= param->rTol * (32.0 / sqrt(3.0))) {
      Fatal("Maximum segment length must be > rTol * 32 / sqrt(3.0)\n"
            "    Current maxSeg = %lf, rTol = %lf",
            param->maxSeg, param->rTol);
    }

    if (param->minSeg > (0.5 * param->maxSeg)) {
      Fatal("Minimum segment length must be < (0.5 * maxSeg)\n"
            "    Current minSeg = %lf, maxSeg = %lf",
            param->minSeg, param->maxSeg);
    }

    if (param->maxSeg <= param->minSeg) {
      Fatal("Max segment length (%e) must be greater than the\n"
            "    minimum segment length (%e)",
            param->maxSeg, param->minSeg);
    }

    if (param->maxSeg > (minCellSize * 0.9)) {
      Fatal("The maxSeg length must be less than the "
            "minimum cell size * 0.9.  Current values:\n"
            "    maxSeg    = %.1f\n    cellWidth = %.1f",
            param->maxSeg, minCellSize);
    }

    if (param->remeshAreaMin > (0.25 * param->remeshAreaMax)) {
      Fatal("remeshAreaMin must be less than 0.25*remeshAreaMax\n"
            "    Current remeshAreaMin = %lf, remeshAreaMax = %lf",
            param->remeshAreaMin, param->remeshAreaMax);
    }

    /*
     *          Now check for conditions that although not fatal, may result
     *          in undesired behaviour, and warn the user.
     */
    if (param->rc < 0.1) {
      fprintf(stderr, "WARNING: specified rc value (%e) will "
                      "yield a \nnegative core energy\n",
              param->rc);
    }

    tmp = (param->maxSeg * param->maxSeg * param->maxSeg);

    if (param->remeshAreaMax > (0.25 * sqrt(3.0) * tmp)) {
      fprintf(stderr, "WARNING: Area criteria will be unused "
                      "in remesh operations!\n");
      fprintf(stderr, "         rmeshAreaMax = %lf, maxSeg = %lf\n",
              param->remeshAreaMax, param->maxSeg);
    }

    if (param->rann > (0.5 * param->rc + eps)) {
      fprintf(stderr, "WARNING: Separation distance is larger "
                      "than the core radius!\n");
      fprintf(stderr, "         rann = %lf, rc = %lf\n", param->rann,
              param->rc);
    }

    if (param->rann > (2.0 * param->rTol)) {
      fprintf(stderr, "WARNING: Collision distance is outside the "
                      "position error tolerance!\n");
      fprintf(stderr, "         rann = %lf, rTol = %lf\n", param->rann,
              param->rTol);
    }
  } /* if domain == 0 */

  /*
   *      Based on the mobility law selected in the control
   *      parameter file, set:
   *        1) the material type (BCC, FCC, etc.)
   *        2) the specific mobility type
   *        3) a pointer to the proper mobility function
   *        4) number of burgers vector groups used in
   *           tracking dislocation density per burgers vector
   *
   *      *************************************************
   *      ***                                           ***
   *      ***                  IMPORTANT!               ***
   *      ***   If you change any numBurgGroups value   ***
   *      ***   specified below, you must change the    ***
   *      ***   DENSITY_FILE_VERSION number defined     ***
   *      ***   in WriteProp.c!                         ***
   *      ***                                           ***
   *      *************************************************
   */
  param->materialType = MAT_TYPE_FCC;
  param->mobilityType = MOB_FCC_0;
  param->mobilityFunc = Mobility_FCC_0;
  param->numBurgGroups = 7;

  /*
   *      If type 1 domain decompositionis enabled, the DLBfreq
   *      value MUST be a multiple of 3.
   */
  if (param->DLBfreq > 0) {
    param->DLBfreq = (param->DLBfreq + 2) / 3 * 3;
  }

  /*
   *      Several portions of the code need to calculate the total
   *      volume of the simulation and a volume factor used for
   *      converting dislocation length to density, so set those
   *      factors now.
   */

  param->simVol =
      param->Dimensions[X] * param->Dimensions[Y] * param->Dimensions[Z];
  param->burgVolFactor =
      1.0 / (param->burgMag * param->burgMag * param->simVol);

#ifdef FULL_N2_FORCES
  /*
   *      To do full n^2 force calculations without remote forces, we need
   *      to explicitly set some flags.
   */
  param->fmEnabled = 0;
  param->DLBfreq = 0;
#endif
}

/*---------------------------------------------------------------------------
 *
 *      Function:     CheckForGlidePlanes
 *      Description:  Some of the mobility modules enforce the motion
 *                    of dislocation segments along specific glide planes.
 *                    If such a mobility is being used, verify that a
 *                    glide plane has been specified for each segment and
 *                    abort with an error if there are any segments without.
 *
 *-------------------------------------------------------------------------*/
static void CheckForGlidePlanes(Home_t *home) {
  int i, j;
  int haveGlidePlanesLocal = 1, haveGlidePlanesGlobal = 0;
  real8 tmp, eps = 1.0e-03;
  Node_t *node;
  Param_t *param;

  param = home->param;

  /*
   *      Loop through every node and every segment attached to the
   *      node and look for zeroed glide planes.
   */
  for (i = 0; i < home->newNodeKeyPtr; i++) {

    if ((node = home->nodeKeys[i]) == (Node_t *)NULL) {
      continue;
    }

    for (j = 0; j < node->numNbrs; j++) {
      tmp = (node->nx[j] * node->nx[j]) + (node->ny[j] * node->ny[j]) +
            (node->nz[j] * node->nz[j]);
      if (tmp < eps) {
        haveGlidePlanesLocal = 0;
        break;
      }
    }

    if (haveGlidePlanesLocal == 0) {
      break;
    }
  }

  /*
   *      Find out if any procesor located a segment with a glide plane
   */

  DSMPI::AllReduce(&haveGlidePlanesLocal, &haveGlidePlanesGlobal, 1, MPI_INT,
                   MPI_MAX);

  /*
   *      If there were any segments with specified glide planes
   *      have processor zero print an error message and force
   *      an abort.
   */
  if ((home->myDomain == 0) && (haveGlidePlanesGlobal == 0)) {
    Fatal("The selected mobility law (%s) requires glide\n"
          "       planes to be defined for all segments.  One or\n"
          "       more segments in the provided nodal data file\n"
          "       do not have glide planes specified.",
          param->mobilityLaw);
  }
}

/*---------------------------------------------------------------------------
 *
 *      Author:       Gregg Hommes
 *
 *      Function:     VerifyBurgersVectors
 *
 *      Description:  This function does a simple sanity check during
 *                    initialization to verify that the burgers vector
 *                    at one end of a segment matches the burgers vector
 *                    at the other end.  Note: This function assumes the
 *                    local domain has nodal info for all nodes terminating
 *                    segments attached to local nodes.  This means that
 *                    the ghost node information must be available before
 *                    this function is called.
 *
 *      Last Modified:  01/09/08: - original version
 *
 *-------------------------------------------------------------------------*/
static void VerifyBurgersVectors(Home_t *home) {
  int nodeID, armID, nbrArmID;
  real8 burgSumX, burgSumY, burgSumZ;
  real8 eps = 1.0e-03;
  Node_t *node, *nbr;

  /*
   *      Loop over all local nodes
   */
  for (nodeID = 0; nodeID < home->newNodeKeyPtr; nodeID++) {

    node = home->nodeKeys[nodeID];

    if (node == (Node_t *)NULL) {
      continue;
    }

    /*
     *          Loop over every segment attached to the local node
     */
    for (armID = 0; armID < node->numNbrs; armID++) {

      nbr = GetNeighborNode(home, node, armID);
      // cout <<"Debug" << armID <<"\t" << node->numNbrs <<endl;
      if (nbr == (Node_t *)NULL) {
        Fatal("VerifyBurgersVectors(): Lookup of node "
              "(%d,%d) failed!",
              node->nbrTag[armID].domainID, node->nbrTag[armID].index);
      }

      /*
       *              Find the neighbor's arm that connects back to the current
       *              node and get its index
       */
      for (nbrArmID = 0; nbrArmID < nbr->numNbrs; nbrArmID++) {
        if ((nbr->nbrTag[nbrArmID].domainID == home->myDomain) &&
            (nbr->nbrTag[nbrArmID].index == node->myTag.index)) {
          break;
        }
      }

      if (nbrArmID >= nbr->numNbrs) {
        cout << node->myTag.domainID << "\t" << node->myTag.index << "\t"
             << nbr->myTag.domainID << "\t" << nbr->myTag.index << endl;
        // cout <<"Coardiantes of node" << node->x <<"\t" <<node->y <<"\t"
        // <<node->z<<endl; cout <<"Coardiantes of neigbor" << nbr->x <<"\t"
        // <<nbr->y <<"\t" <<nbr->z<<endl; cout <<"IDs" <<"\t" << nbrArmID <<
        // "\t" << nbr->numNbrs <<endl;
        Fatal("VerifyBurgersVectors(): neighbor node (%d,%d) "
              "not linked back\n    to local node (%d,%d)",
              nbr->myTag.domainID, nbr->myTag.index, node->myTag.domainID,
              node->myTag.index);
      }

      /*
       *              If the sum of any of the corresponding components of the
       *              burgers vectors at the two ends of the segment are not
       *              equal, we have a problem.
       */
      burgSumX = node->burgX[armID] + nbr->burgX[nbrArmID];
      burgSumY = node->burgY[armID] + nbr->burgY[nbrArmID];
      burgSumZ = node->burgZ[armID] + nbr->burgZ[nbrArmID];

      if ((fabs(burgSumX) > eps) || (fabs(burgSumY) > eps) ||
          (fabs(burgSumZ) > eps)) {
        Fatal("VerifyBurgersVectors(): Burgers vector mismatch!\n"
              "    Segment (%d,%d)--(%d,%d)\n"
              "    burg at first node  = %e %e %e\n"
              "    burg at second node = %e %e %e\n",
              node->myTag.domainID, node->myTag.index, nbr->myTag.domainID,
              nbr->myTag.index, node->burgX[armID], node->burgY[armID],
              node->burgZ[armID], nbr->burgX[nbrArmID], nbr->burgY[nbrArmID],
              nbr->burgZ[nbrArmID]);
      }
    }
  }

  return;
}

/*-------------------------------------------------------------------------
 *
 *      Function:    PrintBanner
 *      Description: Prints a banner to stdout with some basic info
 *                   about the current execution.
 *
 *------------------------------------------------------------------------*/
static void PrintBanner(Home_t *home) {
  int i;
  uid_t uid;
  time_t currTime;
  char workingDir[512];
  char currTimeStr[128];
  char userName[32];
  struct utsname utsInfo;
  struct passwd *pwInfo;

  time(&currTime);
  strcpy(currTimeStr, ctime(&currTime));
  currTimeStr[strlen(currTimeStr) - 1] = 0;

  (void *)getcwd(workingDir, sizeof(workingDir) - 1);
  (void)uname(&utsInfo);

  uid = getuid();

  printf("***********************************************************\n");
  printf("**** \n");
  printf("**** Paradis JHU Running\n");
  printf("**** Modified by Ahmed M. Hussein ... 2013\n");
  printf("**** Time of run:     %s\n", currTimeStr);
  printf("**** Working dir:     %s\n", workingDir);
  printf("**** Execution host:  %s (%s %s %s)\n", utsInfo.nodename,
         utsInfo.sysname, utsInfo.release, utsInfo.machine);
  printf("**** Number of tasks: %d\n", home->numDomains);
  printf("**** \n");
  printf("***********************************************************\n");
  fflush(NULL);
}

/*---------------------------------------------------------------------------
 *
 *      Function:     Initialize
 *      Description:  This is the driver routine for initialization,
 *                    handling some of the general initializations and
 *                    calling all the more specific initialization routines.
 *
 *      Last Modified:  01/09/08: Gregg Hommes - added call to
 *                                VerifyBurgersVectors() as a sanity check.
 *
 *-------------------------------------------------------------------------*/
void Initialize(Home_t *home, const string &sControlFileName, bool bIsFirst) {
  int i;
  char *sep, *start;
  char tmpDataFile[256], tmpStressFile[256];
  char *dataFile, *stressFile;
  Param_t *param;
  InData_t *inData;
  int fd;

  /*
   *      Initialize map between old and new node tags before
   *      reading nodal data and distributing it to the remote
   *      domains.
   */
  home->tagMap = (TagMap_t *)NULL;
  home->tagMapSize = 0;
  home->tagMapEnts = 0;

  home->param = (Param_t *)calloc(1, sizeof(Param_t));
  param = home->param;

  inData = (InData_t *)calloc(1, sizeof(InData_t));

  /*
   *      Verify the command line syntax and pull out the control and
   *      data file names (if provided).  If no control file is specified,
   *      use a default file name If no dataFile name is provided, use
   *      the control file name with the file suffix changed to the
   *      appropriate data file name suffix.  Only need to do this on
   *      processor zero.
   */
  if (home->myDomain == 0) {
    dataFile = NULL;
    // form the data file name
    start = NULL;
    strcpy(tmpDataFile, sControlFileName.c_str());
    start = strrchr(tmpDataFile,
                    '/'); // look for the start of the control file name
    if (start == NULL) {
      start = tmpDataFile;
    }
    sep = strrchr(start, '.');
    // look for the period "." in the control file name

    if ((sep != NULL) && (sep > start)) {
      *sep = 0; // cut the file name starting the period
    }

    strcat(tmpDataFile,
           NODEDATA_FILE_SUFFIX); // concatenate ".data" to the file name
    dataFile = tmpDataFile;       // this becomes the data file name
    // Yejun
    stressFile = NULL;
    start = NULL;
    strcpy(tmpStressFile, sControlFileName.c_str());
    start = strrchr(tmpStressFile,
                    '/'); // look for the start of the control file name
    if (start == NULL) {
      start = tmpStressFile;
    }
    sep = strrchr(start, '.');
    if ((sep != NULL) && (sep > start)) {
      *sep = 0; // cut the file name starting the period
    }
    strcat(tmpStressFile, STRESS_FILE_SUFFIX);
    stressFile = tmpStressFile; // print the banner
    if (bIsFirst) {
      PrintBanner(home);
    }

    // allocate param list memory
    home->ctrlParamList = (ParamList_t *)calloc(1, sizeof(ParamList_t));
    home->dataParamList = (ParamList_t *)calloc(1, sizeof(ParamList_t));
    // --Yejun
    home->stressParamList = (ParamList_t *)calloc(1, sizeof(ParamList_t));

    CtrlParamInit(param, home->ctrlParamList);
    DataParamInit(param, home->dataParamList);
    // Yejun
    StressParamInit(param, home->stressParamList);
    /*
     *          Read runtime parameters from the control file.
     */
    printf("Initialize: Parsing control file %s\n", sControlFileName.c_str());
    ReadControlFile(home, sControlFileName.c_str());
    printf("Initialize: Control file parsing complete\n");

    /*
     *          Some checks on input consistency
     */
    if (bIsFirst) {
      InputSanity(home);
    }
  } /* if (home->myDomain == 0) */

  /*
   *      All domains need the current Param_t structure that's
   *      been populated by domain zero, so broadcast it out.  There
   *      may be subsequent changes to the Param_t structure, but if
   *      so, those updates will be distributed as needed.
   */
  DSMPI::Broadcast((char *)param, sizeof(Param_t), MPI_CHAR, 0);

  /*
   *      Each process needs to get some information about which IO group
   *      it is in before the nodal data files are read in.
   */
  GetParallelIOGroup(home);

  /*
   *      Read the nodal data (and associated parameters) from the
   *      data file (which may consist of multiple file segments).
   */
  ReadNodeDataFile(home, inData, dataFile);
  // update the data file version for writing
  param->dataFileVersion = 6;

  /*
   *         Yejun. Try to open the stress data file.
   */
  if (param->Ttype > 0 || param->SType > 0) {
    ReadStressDataFile(home, inData, stressFile);
    // printf("it is %d %d %lf %lf\n",home->myDomain,
    // home->stress[0][0][0][0].n_x, home->stress[1][1][1][0].x,
    // home->stress[2][2][2][0].s_xy);
  }
  /*
   *      Calculate the length of the problem space in each
   *      of the dimensions
   */
  /*param->Dimensions[0] = param->Dimensions[0]-param->center[0];
      param->Dimensions[1] = param->Dimensions[1]-param->center[1];
      param->Dimensions[2] = param->Dimensions[2]-param->center[2];
      param->Dimensions[3] = param->Dimensions[3]-param->center[0];
      param->Dimensions[4] = param->Dimensions[4]-param->center[1];
      param->Dimensions[5] = param->Dimensions[5]-param->center[2];
      param->Dimensions[0]= param->Dimensions[3]- param->Dimensions[0];
      param->Dimensions[1]= param->Dimensions[4]- param->Dimensions[1];
      param->Dimensions[2]= param->Dimensions[5]- param->Dimensions[2];*/
  param->Lx = param->Dimensions[X];
  param->Ly = param->Dimensions[Y];
  param->Lz = param->Dimensions[Z];
  cout << "C" << param->center[0] << endl;
  cout << "C" << param->center[1] << endl;
  cout << "C" << param->center[2] << endl;
  cout << "Lx" << param->Lx << endl;
  cout << "Ly" << param->Ly << endl;
  cout << "Lz" << param->Lz << endl;
  cout << "D1" << param->Dimensions[0] << endl;
  cout << "D2" << param->Dimensions[1] << endl;
  cout << "D3" << param->Dimensions[2] << endl;
  /*cout << "D4" <<param->Dimensions[3]<<endl;
cout << "D5" <<param->Dimensions[4]<<endl;
cout << "D6" <<param->Dimensions[5]<<endl;*/
  param->invLx = 1.0 / param->Lx;
  param->invLy = 1.0 / param->Ly;
  param->invLz = 1.0 / param->Lz;

  /*
   *      Now that the param structure is fully populated, do any
   *      remaining sanity checks or initializations that could not
   *      or have not been done yet.
   */
  // if(bIsFirst)
  //{
  SetRemainingDefaults(home);
  //}

  /*
   *      Some of the control file parameters are only used when
   *      specific other parameters/toggles have been enabled for
   *      the simulation.  Here we try to identify parameters that
   *      are not used in the current simulation and flag them
   *      so they will not be written into the restart files.  Helps
   *      remove some of the clutter.
   *
   *      This only needs to be done on domain 0, but requires
   *      SetRemainingDefaults() to have be called first to complete
   *      initializations.
   */

  if (home->myDomain == 0) {
    DisableUnneededParams(home);
  }

  /*
   *      Some of the mobility modules require glides planes to be
   *      defined for all segments.  Check for those here.
   */
  CheckForGlidePlanes(home);
  /*
   *      Free up some of the temporary buffers that *may* have been allocated
   */
  FreeInitArrays(home, inData);
  free(inData);

  /*
   *      We attempted to preserve the node tags from the previous
   *      run, so the nodeKeys array may be sparsely populated.  If
   *      that is the case, we have to add all unused tags lower than
   *      <newNodeKeyPtr> to the recycled node array or those tags
   *      will never be used again.
   */
  InitRecycleNodeHeap(home);

  /*
   *      Find out which cells intersect this domain (the native cells), and
   *      their neighbors cells. Find out which domains intersect each of these
   *      native and ghost cells.
   */
  InitCellNatives(home);
  InitCellNeighbors(home);
  InitCellDomains(home);

  /*
   *      For each neighbor domain, build a list of native cells to send to that
   *      domain.
   */
  InitRemoteDomains(home);

  /*
   *      Each domain still needs to communicate with its neighboring domains
   *      to map old arm tags to new ones for any nodes that were retagged
   *      during initialization.
   */
  DistributeTagMaps(home);

/*
 *      Allocate an array to store the cell charge tensor for each cell
 *      (but don't do it if remote forces are being disabled by the
 *      FULL_N2_FORCES flag)
 */
#ifndef FULL_N2_FORCES
  home->cellCharge = (real8 *)malloc(param->nXcells * param->nYcells *
                                     param->nZcells * 9 * sizeof(real8));
#endif

  /*
   *      Initialize operation list used for collisions and other topological
   *      changes.
   */
  InitOpList(home);

#ifndef FULL_N2_FORCES
/*
 *      If the Fast Multipole code is enabled, initialize the image
 *      correction table.  Otherwise (if necessary) read in PBC image
 *      correction table and PBC stress tables BEFORE creation of and
 *      cd to the output file directory.
 *
 *      NOTE: This call to CorrectionTableInit() MUST be done after
 *      the first call to FMInit() (which is invoked from within the
 *      InitCellNatives() function called above).
 */
// 		ReadRijm(home);
// 		ReadRijmPBC(home);
#endif

  /*
   *      Create the output directory (and sub-directories) and reset the
   *      current working directory to the top level output directory.
   */
  if (home->param->dirname[0] != 0) {
    OpenDir(home);
  }

  /*
   *      Do the initial sort of native nodes into their proper subcells
   *      and send the initial ghost node data.  Previously this was at
   *      the top of the ParadisStep loop, but modifications to the way
   *      output (including X-windows plot data) is generated requires that
   *      the cpus have accurate ghost node data before calling GenerateOutput.
   *      Hence, we do the first ghost node communications here, and move
   *      the subsequent ghost node comm from the beginning of the ParadisStep
   *      loop to the end.
   */
  SortNativeNodes(home);
  CommSendGhosts(home);

  /*
   *      Have each processor look at all segments attached to its local
   *      nodes and verify that the burgers vector is the same at both
   *      ends of the node.  Just a sanity check to prevent people from
   *      doing silly things like creating a nodal configuration by hand
   *      and putting in inconsistent burgers vectors.
   */
  VerifyBurgersVectors(home);

  /*
   *      If necessary, create the rotation matrices needed for rotating
   *      from geometries defined in the users laboratory frame to the
   *      standard crystalographic frame and back.
   */
  if (param->useLabFrame) {

    Normalize(&param->labFrameXDir[0], &param->labFrameXDir[1],
              &param->labFrameXDir[2]);

    Normalize(&param->labFrameYDir[0], &param->labFrameYDir[1],
              &param->labFrameYDir[2]);

    NormalizedCrossVector(param->labFrameXDir, param->labFrameYDir,
                          param->labFrameZDir);

    for (i = 0; i < 3; i++) {
      home->rotMatrix[0][i] = param->labFrameXDir[i];
      home->rotMatrix[1][i] = param->labFrameYDir[i];
      home->rotMatrix[2][i] = param->labFrameZDir[i];
    }

    Matrix33Invert(home->rotMatrix, home->rotMatrixInverse);
  }

  CheckMemUsage(home, "Initialize");
}
