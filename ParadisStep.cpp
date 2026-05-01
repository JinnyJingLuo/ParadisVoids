/*-------------------------------------------------------------------------
 *
 *      Function:     ParadisStep
 *      Description:  This function controls everything needed for a
 *                    single step of a ParaDiS simulation including
 *                    force calculations, ghost cell communications,
 *                    node migration, dynamic load balance, output
 *                    generation, etc.
 *
 *-----------------------------------------------------------------------*/

#include <stdio.h>
#include <iostream> //added
#include <fstream>  //added
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include "Home.h"
#include "Util.h"
#include "Comm.h"
#include "Mobility.h"
#include "Decomp.h"
#include "QueueOps.h"
#include "ParadisCrossSlipServer.h"
#include "DSMPI.h"
#include "time.h"
#include <ctime>
#include "TwinPlaneCrossSlip.h" //junjie
//#include "Windows.h" //added
using namespace std; // added
                     // void timer();

/*
 *      By default, there are no runtime checks to see if all of
 *      the dislocations have annihilated themselves.  To enable
 *      a check with an abort if it happens, simply define the
 *      DEBUG_CHECK_FOR_ZERO_SEG value below to 1 rather than zero.
 */
#define DEBUG_CHECK_FOR_ZERO_SEG 0

/*
 *      For debugging only.  If DEBUG_STEP is not defined, all
 *      calls to Synchronize() will be replaced with an empty
 *      block of code, but if it is defined, the calls will
 *      be replaced with a call to syncronize the code and log
 *      a message.
 */
#ifdef DEBUG_STEP
#define Synchronize(a, b) _Synchronize((a), (b))
#else
#define Synchronize(a, b)                                                      \
  {}
#endif

/*
 *      Explicitly synchronize parallel execution via an MPI
 *      barrier and print a message when all tasks have reached
 *      the barrier.  For debug only.
 */
void _Synchronize(Home_t *home, char *msg) {
  // ofstream outFile;
  // outFile.open("synchronize.data");

  // int start_s1=clock();

  DSMPI::Barrier();

  if (home->myDomain == 0) {
    printf(" *** %s: All tasks synchronized\n", msg);
    fflush(NULL);
  }
  // int stop_s1=clock();
  // outFile << "time: " << (stop_s1-start_s1)/double(CLOCKS_PER_SEC)*1000 <<
  // endl; outFile.close(); timer ();//added
}

void PrintTime(Home_t *home, const unsigned int &iStampID) {
  // ofstream outFile;
  // outFile.open("print_time.data");

  // int start_s2=clock();
  if (home->myDomain == 0) {
    printf("time %d : %25.20f\n", iStampID, MPI_Wtime());
    fflush(NULL);
  }
  DSMPI::Barrier();
  // int stop_s2=clock();
  // outFile << "time: " << (stop_s2-start_s2)/double(CLOCKS_PER_SEC)*1000 <<
  // endl; outFile.close(); timer (); //added
}
void WriteVelocity(Home_t *poHome) {
  // ofstream outFile;
  // outFile.open("write_velocity.data");

  // int start_s3=clock();

  // compute the mean velocity
  unsigned int i = 0;
  Node_t *poNode = NULL;
  double dTemp = 0.0;
  double dLocalSums[3] = {0.0, 0.0, 0.0};
  for (i = 0; i < poHome->newNodeKeyPtr; i++) {
    poNode = poHome->nodeKeys[i];
    if (poNode == NULL) {
      continue;
    }
    dLocalSums[0] = dLocalSums[0] + 1.0; // the total number of nodes
    dTemp = poNode->vX * poNode->vX + poNode->vY * poNode->vY +
            poNode->vZ * poNode->vZ;
    dLocalSums[2] = dLocalSums[2] + dTemp; // the sum of squared velocities
    dTemp = sqrt(dTemp);
    dLocalSums[1] = dLocalSums[1] + dTemp; // the sum of velocities
  }
  double dGlobalSums[3] = {0.0, 0.0, 0.0};
  DSMPI::AllReduce(dLocalSums, dGlobalSums, 3, MPI_DOUBLE, MPI_SUM);
  unsigned int iNodesCount = (int)floor(dGlobalSums[0] + 0.5);
  double dVelocityMean = dGlobalSums[1] / dGlobalSums[0];
  double dVelocityVariance =
      dGlobalSums[2] / dGlobalSums[0] - dVelocityMean * dVelocityMean;
  FILE *fpFile = NULL;
  if ((poHome->myDomain == 0) &&
      (poHome->cycle % poHome->param->savefreq == 0)) {
    fpFile = fopen("velocity_stats", "a");
    fprintf(fpFile, "%d : %d\t\t%e\t\t%e\n", poHome->cycle, iNodesCount,
            dVelocityMean, sqrt(dVelocityVariance));
    fclose(fpFile);
  }
  DSMPI::Barrier();
  // timer (); //added
  // int stop_s3=clock();
  // outFile << "time: " << (stop_s3-start_s3)/double(CLOCKS_PER_SEC)*1000 <<
  // endl; outFile.close();
}

void InitializeStep(Home_t *poHome) {
  // ofstream outFile;
  // outFile.open("Initialize_time.data");

  // int start_s4=clock();
  // increment the simulation cycle
  poHome->cycle++;
  // increment the time step
  poHome->param->timeNow = poHome->param->timeNow + poHome->param->deltaTT;
  // Yejun
  if (poHome->param->Ttype > 0 || poHome->param->SType > 0) {
    if (poHome->stress[0][0][0][poHome->param->Time_evo - 1].t_end <
            poHome->param->timeNow ||
        poHome->stress[0][0][0][0].t_start > poHome->param->timeNow) {
      poHome->param->stress_timestep = -1;
    } else {
      for (int i = 0; i < poHome->param->Time_evo; i++) {
        if ((poHome->stress[0][0][0][i].t_start <= poHome->param->timeNow) &&
            (poHome->stress[0][0][0][i].t_end > poHome->param->timeNow)) {
          poHome->param->stress_timestep = i + 1;
          // printf("%d", poHome->param->stress_timestep);
          // poHome->param->TempK = poHome->stress[0][0][0][i].Temp;
          break;
        }
      }
    }
  }
  // print current cycle information
  if (poHome->myDomain == 0) {
    time_t tp;
    time(&tp);
    if ((poHome->param->Ttype > 0 || poHome->param->SType > 0) &&
        poHome->param->stress_timestep > 0) {
      printf(
          "cycle=%-8d  timestep=%e  timeNow=%e  stress_time=%e  %s",
          poHome->cycle, poHome->param->deltaTT, poHome->param->timeNow,
          poHome->stress[0][0][0][poHome->param->stress_timestep - 1].t_start,
          asctime(localtime(&tp)));
    } else {
      printf("cycle=%-8d  timestep=%e  timeNow=%e  %s", poHome->cycle,
             poHome->param->deltaTT, poHome->param->timeNow,
             asctime(localtime(&tp)));
    }
  }
  // increment the save counter if this is a writing step
  if (poHome->cycle % poHome->param->savefreq == 0) {
    poHome->param->savecounter++;
    sprintf(poHome->param->APBFileName, "apb_%04d.dat",
            poHome->param->savecounter);
  }
  poHome->param->mergedelpStrain[0] = 0.0;
  poHome->param->mergedelpStrain[1] = 0.0;
  poHome->param->mergedelpStrain[2] = 0.0;
  poHome->param->mergedelpStrain[3] = 0.0;
  poHome->param->mergedelpStrain[4] = 0.0;
  poHome->param->mergedelpStrain[5] = 0.0;

  poHome->param->mergedelpSpin[0] = 0.0;
  poHome->param->mergedelpSpin[1] = 0.0;
  poHome->param->mergedelpSpin[2] = 0.0;
  poHome->param->mergedelpSpin[3] = 0.0;
  poHome->param->mergedelpSpin[4] = 0.0;
  poHome->param->mergedelpSpin[5] = 0.0;
  // timer (); //added
  // int stop_s4=clock();
  // outFile << "time: " << (stop_s4-start_s4)/double(CLOCKS_PER_SEC)*1000 <<
  // endl; outFile.close();
}

void ReportForce(Home_t *home) {
  // ofstream outFile;
  // outFile.open("Report_force.data");

  // int start_s5=clock();
  unsigned int i = 0;
  Node_t *poNode = NULL;
  for (i = 0; i < home->newNodeKeyPtr; i++) {
    poNode = home->nodeKeys[i];
    if (poNode == NULL)
      continue;
    if (poNode->constraint >=SURFACE_NODE) {
      printf("force on surface node (%d,%d) : %e,%e,%e\n",
             poNode->myTag.domainID, poNode->myTag.index, poNode->fX,
             poNode->fY, poNode->fZ);
      printf("surface normal of surface node : %e,%e,%e\n", poNode->dNx,
             poNode->dNy, poNode->dNz);
    }
  }
  // timer ();//added
  // int stop_s5=clock();
  // outFile << "time: " << (stop_s5-start_s5)/double(CLOCKS_PER_SEC)*1000 <<
  // endl; outFile.close();
}

void ParadisStep(Home_t *home) {
  // ofstream outFile;
  // outFile.open("net_charge_tensor.data");

  // int start_s6=clock();
  int i;
  Param_t *param;
  Node_t *node;
  param = home->param;
  InitializeStep(home);

  /*
   *      Calculate the net charge tensor for each cell (includes global comm)
   */
  CellCharge(home);
  // int stop_s6=clock();
  // outFile << "time: " << (stop_s6-start_s6)/double(CLOCKS_PER_SEC)*1000 <<
  // endl; outFile.close();

  /*
   *      Define load curve and calculate the increment in applied stress this
   * cycle
   */
  // outFile.open("load_curve.data");
  // int start_s7=clock();
  LoadCurve(home);
  // int stop_s7=clock();
  // outFile << "time: " << (stop_s7-start_s7)/double(CLOCKS_PER_SEC)*1000 <<
  // endl; outFile.close();

  /*
   *      Calculate new force and velocity data for all nodes or a selected
   *      subset and distribute the new data out to neighboring domains.
   */

  // outFile.open("new_force_velocity.data");
  // int start_s8=clock();
  NodeForce(home, home->poExternalLoadServer, home->poPrecipitateServer);
  CalcNodeVelocities(home);
  CommSendVelocity(home);
  WriteVelocity(home);
  // int stop_s8=clock();
  // outFile << "time: " << (stop_s8-start_s8)/double(CLOCKS_PER_SEC)*1000 <<
  // endl; outFile.close();
  /*
   *      Invoke the selected time step integration method.  The
   *      selected method will calculate the time step as well as
   *      move nodes to their correct locations and communicate
   *      the new nodal force/velocity data to neighboring domains.
   */

  // outFile.open("Euler_Integartor.data");
  // int start_s9=clock();
  ForwardEulerIntegrator(home);
  // handle the precipitates, if any
  if (home->param->ShearablePrecipitates == 0) {
    home->poPrecipitateServer->CheckNodes(home);
  } else {
    home->poPrecipitateServer->UpdateAPB(home);
  }
  DSMPI::Barrier();
  if (home->param->EnableTwinPlaneCrossSlip == 1) {
    // junjie
    ClearOpList(home);
    TwinPlaneCrossSlip(home);
    DSMPI::Barrier();
    CommSendRemesh(home);
    FixRemesh(home);
    // junjie
  }

  if (home->param->BoundaryType != PERIODIC_BOUNDARY) {
    home->poSurface->CheckNodes(home);
    DSMPI::Barrier();
    home->poSurface->StoreSurfaceArms(home);
    home->poSurface->StoreSurfaceNodesMotion(home);
    if (home->cycle % home->param->savefreq == 0) {
      home->poSurface->WriteSurfaceSegments(home);
    }
  }
  DSMPI::Barrier();

  // junjie
  // int stop_s9=clock();
  // outFile << "time: " << (stop_s9-start_s9)/double(CLOCKS_PER_SEC)*1000 <<
  // endl; outFile.close();
  /*
   *      Increment the per-burgers vector density gain/loss with
   *      changes for this cycle.  This must be done immediately
   *      after timestep integration!
   *
   *      Note: This is currently only applicable to BCC simulations.
   */

  // outFile.open("getDensityDElta.data");

  // int start_s10=clock();

  GetDensityDelta(home);
  // int stop_s10=clock();
  // outFile << "time: " << (stop_s10-start_s10)/double(CLOCKS_PER_SEC)*1000 <<
  // endl; outFile.close();
  /*
   *      Calculate the new plastic strain.
   */
  // outFile.open("new_plasticStrain.data");

  // int start_s11=clock();
  DeltaPlasticStrain(home);
  ParadisCrossSlipServer::GetInstance()->HandleCrossSlip(home);
  if (home->param->BoundaryType != PERIODIC_BOUNDARY &&
      home->param->EnableTwinPlaneCrossSlip != 1) {
    home->poSurface->CheckNodes(home);
    DSMPI::Barrier();
  }
  // int stop_s11=clock();
  // outFile << "time: " << (stop_s11-start_s11)/double(CLOCKS_PER_SEC)*1000 <<
  // endl; outFile.close();
  /*
   *      The call to GenerateOutput will update the time and cycle counters,
   *      determine if any output needs to be generated at this stage, and
   *      call the appropriate I/O functions if necessary.
   */

  // outFile.open("GenerateOutput.data");

  // int start_s12=clock();

  GenerateOutput(home, STAGE_CYCLE);
  // int stop_s12=clock();
  // outFile << "time: " << (stop_s12-start_s12)/double(CLOCKS_PER_SEC)*1000 <<
  // endl; outFile.close();
  /*
   *      Before doing topological changes, set flags indicating any
   *      nodes exempt from topological changes.  These flags are used
   *      in both splitting multi-arm nodes and collisions, so this
   *      function should be invoked before either of those items are done.
   */
  // outFile.open("topologFlag.data");

  // int start_s13=clock();

  InitTopologyExemptions(home);
  // junjie
  ExemptCollisionAfterDissociation(home);
  // junjie
  // int stop_s13=clock();
  // outFile << "time: " << (stop_s13-start_s13)/double(CLOCKS_PER_SEC)*1000 <<
  // endl; outFile.close();

  /*
   *      Now do all the topological changes from segment interactions
   *      (collisions, multinode splitting)...  Clear the list of local
   *      operations that will be sent to the remote domains for processsing,
   *      then split any multi-arm nodes that need splitting, cross slip
   *      nodes (as needed/allowed), handle all local collisions, then
   *      send remote nodes the list of ops needed to keep their data in sync.
   */

  // outFile.open("topologicalChanges.data");

  // int start_s14=clock();
  ClearOpList(home);
  SortNodesForCollision(home);
  // int stop_s14=clock();
  // outFile << "time: " << (stop_s14-start_s14)/double(CLOCKS_PER_SEC)*1000 <<
  // endl; outFile.close();

  /*
   *      Search for dislocation segments in close proximity to each other
   *      and if necessary handle any collision between them.
   */

  // outFile.open("Dislocation_Proximity.data");

  // int start_s15=clock();
  home->poCollisionServer->HandleCollisions(home);

  DSMPI::Barrier();

  CommSendRemesh(home);
  FixRemesh(home);
  // int stop_s15=clock();
  // outFile << "time: " << (stop_s15-start_s15)/double(CLOCKS_PER_SEC)*1000 <<
  // endl; outFile.close();
  /*
   *      Under certain circumstances, parallel topological changes can
   *      create double links between nodes; links which can not be detected
   *      until after FixRemesh() is called... so, a quick check has to be
   *      done to clean up these potential double-links here, or they will
   *      cause problems later on.  Should only have to check nodes local
   *      to this domain.
   */
  // outFile.open("parallel_toplogicalChange.data");

  // int start_s16=clock();
  for (i = 0; i < home->newNodeKeyPtr; i++) {
    if ((node = home->nodeKeys[i]) == NULL)
      continue;
    (void)RemoveDoubleLinks(home, node, 0);
    node->flags &= ~NODE_CHK_DBL_LINK;
  }
  // int stop_s16=clock();
  // outFile << "time: " << (stop_s16-start_s16)/double(CLOCKS_PER_SEC)*1000 <<
  // endl; outFile.close();
  if (home->param->EnableTwinPlaneCrossSlip == 1) {
    // junjie
    ClearOpList(home);
    DissociateFourArmsNode(home);
    DSMPI::Barrier();
    CommSendRemesh(home);
    FixRemesh(home);
    // junjie
  }

  /*
   *      Invoke mesh coarsen/refine
   */
  // outFile.open("remesh.data");

  // int start_s17=clock();
  Remesh(home);
  // int stop_s17=clock();
  // outFile << "time: " << (stop_s17-start_s17)/double(CLOCKS_PER_SEC)*1000 <<
  // endl; outFile.close();
  /*
   *      If necessary, use the current load data to generate a new
   *      domain decomposition to rebalance the workload among the
   *      processors.
   */
  // outFile.open("domain_decomposition.data");

  // int start_s18=clock();
  Rebalance(home);
  // int stop_s18=clock();
  // outFile << "time: " << (stop_s18-start_s18)/double(CLOCKS_PER_SEC)*1000 <<
  // endl; outFile.close();
  /*
   *      Send any nodes that have moved beyond the domain's
   *      boundaries to the domain the node now belongs to.
   */
  // outFile.open("Migrate.data");

  // int start_s19=clock();
  Migrate(home);
  // int stop_s19=clock();
  // outFile << "time: " << (stop_s19-start_s19)/double(CLOCKS_PER_SEC)*1000 <<
  // endl; outFile.close();
  /*
   *      Recycle all the ghost nodes: move them back to the free Queue
   */
  // outFile.open("Free_Queue.data");

  // int start_s20=clock();
  RecycleGhostNodes(home);
  // int stop_s20=clock();
  // outFile << "time: " << (stop_s20-start_s20)/double(CLOCKS_PER_SEC)*1000 <<
  // endl; outFile.close();
  /*
   *      Sort the native nodes into their proper subcells.
   */
  // outFile.open("Sorting_nodes.data");

  // int start_s21=clock();
  SortNativeNodes(home);
  // int stop_s21=clock();
  // outFile << "time: " << (stop_s21-start_s21)/double(CLOCKS_PER_SEC)*1000 <<
  // endl; outFile.close();
  /*
   *      Communicate ghost cells to/from neighbor domains
   */
  // outFile.open("Communicate_cells.data");

  // int start_s22=clock();
  CommSendGhosts(home);
  CheckMemUsage(home, "ParadisStep-complete");
  /*int stop_s22=clock();
  outFile << "time: " << (stop_s22-start_s22)/double(CLOCKS_PER_SEC)*1000 <<
  endl; outFile.close();*/
  /*
   *      Zero out the count of force calculations done this cycle
   *      so the load-balancing in the next step is based on accurate
   *      values.
   */
  // outFile.open("count_force.data");

  // int start_s23=clock();
  home->cycleForceCalcCount = 0;
  /*int stop_s23=clock();
  outFile << "time: " << (stop_s23-start_s23)/double(CLOCKS_PER_SEC)*1000 <<
  endl; outFile.close();*/
}
/*
void timer(){
        int sec;
                sec = 0;
                int min;
                min = 0;
                int hour;
                hour = 0;
                while (true){
                        sleep (1000);
                        sec = sec +1;
                        if (sec ==60){
                        min = min+1;
                        sec = 0;
                }
                        if (min==60){
                                hour = hour +1;
                                min =0;
                        }
                        system ("CLS");
                        cout << hour << " hours " << min << " minutes " << " and
" <<sec << " seconds " << endl;
                }

        }




*/
