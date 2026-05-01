/****************************************************************************
 *
 *  Param.h  Define the Param struct, which holds or points to all control
 *           parameters for this simulation
 *
 ***************************************************************************/
#ifndef _PARAM_H
#define _PARAM_H

#include "Parse.h"
#include "Home.h"
#include "Mobility.h"

/*
 *      Define a couple strings related to the nodal data files
 *      written along with the control file containing the global
 *      parameter values.
 */
#define HDF_DATA_FILE_SUFFIX ".hdf"
#define NODEDATA_FILE_SUFFIX ".data"
#define NODEDATA_FILE_VERSION 4
// Suffix for thermal stress data file
#define STRESS_FILE_SUFFIX ".str"
#define PERIODIC_BOUNDARY 0
#define FREE_BOUNDARY 1
#define RIGID_BOUNDARY 2
#define ONLY_PZ_RIGID_BOUNDARY 3

struct _param {
  /*
   *      Defines the number of domains in each dimension of the
   *      entire problem space.
   */
  int nXdoms;
  int nYdoms;
  int nZdoms;

  /*
   *      Defines number of cells in each of the dimensions of the
   *      entire problem space.
   */
  int nXcells;
  int nYcells;
  int nZcells;

  int BoundaryType;

  /*
   *      "natural" min and max cell indices for this domain. (i.e. indices
   *      before being incremented by 1 to allow for ghost cells)
   */
  int iCellNatMin;
  int iCellNatMax;
  int jCellNatMin;
  int jCellNatMax;
  int kCellNatMin;
  int kCellNatMax;

  /*
   *      Domain decomposition and rebalance values
   */
  int DLBfreq; /* how often to load balance */

  int Ttype; /* if the thermal stress data file is used */

  int SType; /* if the thermal stress data file is used */
             /*
              *      Turn on the Peierls stress(thredhold for dislocation glides)
              */
  unsigned int EnablePeierls;
  real8 PeierlsScrew;
  real8 PeierlsEdge;
  real8 PeierlsOthers;

  /*
   *     Turn on twin plane cross slip
   */
  unsigned int EnableTwinPlaneCrossSlip;
  real8 A;
  real8 B;
  real8 C;
  real8 X_1;
  real8 Y_1;
  real8 Z_1; // define the twin plane

  /*
   *      Simulation time and timestepping controls
   */
  int cycleStart; /* Starting cycle number for the simulation */
  int maxstep;    /* Cycles to execute before terminating */
  real8 timeNow;  /* current simulation time */

  real8 deltaTT; /* duration of previous timestep */
  real8 rTol;    /* Maximum error allowed in timestep */
  real8 rmax;    /* maximum migration distance per timestep */
                 /* for any node   */

  /*
   *      Discretization parameters and controls for topological changes
   */
  real8 minSeg; /* min allowable segment length, before */
                /* removing a node */
  real8 maxSeg; /* max allowable segment length, before*/
                /* adding a node*/
  int collisionMethod;
  real8 remeshAreaMax;    /* This is calculated from remeshAreaMin and  */
                          /* remeshAreaRatio and hence not user-provided*/
  real8 remeshAreaMin;    /* This values is based on the minSeg value  */
                          /* and hence not specified by the user.      */
  int splitMultiNodeFreq; /* Code will attempt to split multi-nodes */
                          /* every cycle that is a multiple of this */
                          /* value. */

  /*
   *      Names of tables for non-FMM far-field forces
   */
  char Rijmfile[MAX_STRING_LEN];
  char RijmPBCfile[MAX_STRING_LEN];

  /*
   *      Loading condition parameters
   */
  real8 TempK;            /* Temperature in deg K */
  int loadType;           /* 0 Creep test */
                          /* 1 Constant strain rate test */
                          /* 2 Displacement-controlled test */
                          /* 3 Load-controlled, load vs. time curve */
                          /* 4 Cyclic loading condition */
  real8 appliedStress[6]; /* External stress in units of Pa  */
                          /* as [sigma11, sigma22, sigma33,  */
                          /* sigma23, sigma31, sigma12] when */
                          /*  <loadType> == 0.               */
  real8 eRate;            /* Strain rate. Used when loadType == 1 */
  int indxErate;          /* to be compatible with micro3d */
  real8 edotdir[3];       /* Uniaxial loading direction accompanying */
                          /* eRate                                   */
  real8 cTimeOld;         /* Timestep related to cyclic loading */
  real8 netCyclicStrain;  /* Net accumulated strain under cyclic load */
  real8 dCyclicStrain;    /* Incremental strain under cyclic load */
  int numLoadCycle;       /* Number of cyclic cycles */
  real8 eAmp;             /* Strain amplitude used with cyclic loading */

  // Cross Slip Parameters
  // Added by Ahmed M. Hussein, July 27th 2012
  double BulkCrossSlipActivationEnergy;
  double SurfaceCrossSlipActivationEnergy;
  double HirthCrossSlipActivationEnergy;
  double LCLockCrossSlipActivationEnergy;
  double GlideLockCrossSlipActivationEnergy;

  double BulkCrossSlipActivationVolumeFactor;
  double SurfaceCrossSlipActivationVolumeFactor;
  double HirthCrossSlipActivationVolumeFactor;
  double LCLockCrossSlipActivationVolumeFactor;
  double GlideLockCrossSlipActivationVolumeFactor;

  double BulkCrossSlipReferenceLength;
  double SurfaceCrossSlipReferenceLength;
  double HirthCrossSlipReferenceLength;
  double LCLockCrossSlipReferenceLength;
  double GlideLockCrossSlipReferenceLength;

  double SurfaceCrossSlipLength;
  double HirthCrossSlipLength;
  double LCLockCrossSlipLength;
  double GlideLockCrossSlipLength;

  double BulkCrossSlipFrequency;
  double SurfaceCrossSlipFrequency;
  double HirthCrossSlipFrequency;
  double LCLockCrossSlipFrequency;
  double GlideLockCrossSlipFrequency;

  double BulkStackingFaultEnergy;
  double SurfaceStackingFaultEnergy;
  double HirthStackingFaultEnergy;
  double LCLockStackingFaultEnergy;
  double GlideLockStackingFaultEnergy;

  double BulkDoverBRatio;
  double SurfaceDoverBRatio;
  double HirthDoverBRatio;
  double LCLockDoverBRatio;
  double GlideLockDoverBRatio;

  unsigned int HandleCrossSlip;
  // unsigned int Handle;
  unsigned int CrossSlipHandlingFrequency;

  double CyclicLoadAmplitude;
  double CyclicLoadFrequency;
  double InitialTemperature;
  double HeatingRate;
  double YoungsModulus;
  double PoissonsRatio;
  double MassDensity;
  double ThermalConductivity;
  double ThermalExpansionCoefficient;
  double SpecificHeatCapacity;
  double ReferenceTemperature;
  double FEMTimeStep;

  unsigned int EnableFEM;
  unsigned int FEMCorrectionFrequency;
  unsigned int LoopExtractionFrequency;
  char FEMRestartFileName[MAX_STRING_LEN];

  unsigned int AnnihilateSegments;
  double AnnihilationDistance;

  double APBEnrgy;
  double FrictionStress;
  double SolidSolutionStrength;
  double KWStress;
  double voidR;
  unsigned int VelocityDampingSteps;

  // cross slip statistics
  unsigned int BulkCrossSlipEventsCount;
  unsigned int SurfaceCrossSlipEventsCount;
  unsigned int AttractiveCrossSlipEventsCount;
  unsigned int RepulsiveCrossSlipEventsCount;
  // per slip plane
  unsigned int HandleCrossSlipPerSlipSystem;
  unsigned int SlipSystem1BulkCrossSlipEventsCount;
  unsigned int SlipSystem2BulkCrossSlipEventsCount;
  unsigned int SlipSystem3BulkCrossSlipEventsCount;
  unsigned int SlipSystem4BulkCrossSlipEventsCount;
  unsigned int SlipSystem5BulkCrossSlipEventsCount;
  unsigned int SlipSystem6BulkCrossSlipEventsCount;
  unsigned int SlipSystem7BulkCrossSlipEventsCount;
  unsigned int SlipSystem8BulkCrossSlipEventsCount;
  unsigned int SlipSystem9BulkCrossSlipEventsCount;
  unsigned int SlipSystem10BulkCrossSlipEventsCount;
  unsigned int SlipSystem11BulkCrossSlipEventsCount;
  unsigned int SlipSystem12BulkCrossSlipEventsCount;
  unsigned int OtherSlipSystemBulkCrossSlipEventsCount;
  unsigned int SlipSystem1SurfaceCrossSlipEventsCount;
  unsigned int SlipSystem2SurfaceCrossSlipEventsCount;
  unsigned int SlipSystem3SurfaceCrossSlipEventsCount;
  unsigned int SlipSystem4SurfaceCrossSlipEventsCount;
  unsigned int SlipSystem5SurfaceCrossSlipEventsCount;
  unsigned int SlipSystem6SurfaceCrossSlipEventsCount;
  unsigned int SlipSystem7SurfaceCrossSlipEventsCount;
  unsigned int SlipSystem8SurfaceCrossSlipEventsCount;
  unsigned int SlipSystem9SurfaceCrossSlipEventsCount;
  unsigned int SlipSystem10SurfaceCrossSlipEventsCount;
  unsigned int SlipSystem11SurfaceCrossSlipEventsCount;
  unsigned int SlipSystem12SurfaceCrossSlipEventsCount;
  unsigned int OtherSlipSystemSurfaceCrossSlipEventsCount;
  unsigned int SlipSystem1RepulsiveCrossSlipEventsCount;
  unsigned int SlipSystem2RepulsiveCrossSlipEventsCount;
  unsigned int SlipSystem3RepulsiveCrossSlipEventsCount;
  unsigned int SlipSystem4RepulsiveCrossSlipEventsCount;
  unsigned int SlipSystem5RepulsiveCrossSlipEventsCount;
  unsigned int SlipSystem6RepulsiveCrossSlipEventsCount;
  unsigned int SlipSystem7RepulsiveCrossSlipEventsCount;
  unsigned int SlipSystem8RepulsiveCrossSlipEventsCount;
  unsigned int SlipSystem9RepulsiveCrossSlipEventsCount;
  unsigned int SlipSystem10RepulsiveCrossSlipEventsCount;
  unsigned int SlipSystem11RepulsiveCrossSlipEventsCount;
  unsigned int SlipSystem12RepulsiveCrossSlipEventsCount;
  unsigned int OtherSlipSystemRepulsiveCrossSlipEventsCount;
  unsigned int SlipSystem1AttractiveCrossSlipEventsCount;
  unsigned int SlipSystem2AttractiveCrossSlipEventsCount;
  unsigned int SlipSystem3AttractiveCrossSlipEventsCount;
  unsigned int SlipSystem4AttractiveCrossSlipEventsCount;
  unsigned int SlipSystem5AttractiveCrossSlipEventsCount;
  unsigned int SlipSystem6AttractiveCrossSlipEventsCount;
  unsigned int SlipSystem7AttractiveCrossSlipEventsCount;
  unsigned int SlipSystem8AttractiveCrossSlipEventsCount;
  unsigned int SlipSystem9AttractiveCrossSlipEventsCount;
  unsigned int SlipSystem10AttractiveCrossSlipEventsCount;
  unsigned int SlipSystem11AttractiveCrossSlipEventsCount;
  unsigned int SlipSystem12AttractiveCrossSlipEventsCount;
  unsigned int OtherSlipSystemAttractiveCrossSlipEventsCount;
  double TotalBulkCrossSlippedChainsLength;
  double TotalSurfaceCrossSlippedChainsLength;
  double TotalAttractiveCrossSlippedChainsLength;
  double TotalRepulsiveCrossSlippedChainsLength;

  /*
   *      Values for specifying axes for a user-defined laboratory frame
   */
  int useLabFrame; /* 0 if standard crystalographic frame is to */
                   /* be used, 1 if user-supplied laboratory    */
                   /* frame is used                             */

  real8 labFrameXDir[3]; /* The Z direction is informational only */
  real8 labFrameYDir[3]; /* and is recalculated explicitly from   */
  real8 labFrameZDir[3]; /* the X and Y directions                */

  char mobilityLaw[MAX_STRING_LEN];
  int mobilityType; /* Integer value corresponding to the */
                    /* specified mobility law.  Redundant */
                    /* info, but easier to use in the code*/

  int materialType; /* Type of crystal structure (i.e. BCC, FCC)   */
                    /* This value is set within the code base on   */
                    /* the selected mobility law and is not user-  */
                    /* supplied                                    */

  real8 vacancyConc; /* Concentration of vacancies in the */
                     /* crystal: for calculating osmotic  */
                     /* force (units ?)                   */

  real8 vacancyConcEquilibrium; /* Thermal equilibrium vacacy concen-*/
                                /* ration in defect free crystal: for*/
                                /* calculating osmotic force         */

  real8 shearModulus;
  real8 pois;
  real8 burgMag;

  real8 rc; /* core radius in elastic interaction calculation */

  real8 Ecore; /* core energy (wrt. the choice of rc) in unit of Pa */

  int (*mobilityFunc)(Home_t *home, Node_t *node); /* Set during */
  /* initialization to point to the      */
  /* appropriate mobility function       */

  real8 MobScrew;
  real8 MobEdge;
  real8 MobClimb;

  real8 MobM; /*Yejun: Drag coefficient as a linear function of temp: B=MobM *
                 temp + MobB0 */
  real8 MobB0;

  real8 MobGlide; /* floor on mobility for glide dislocations */
  real8 MobLine;

  /*
   *      Allow for specifying that dislocations with certain types
   *      of burgers vectors or line directions are immobile.
   */
  real8 sessileburgspec[30];
  real8 sessilelinespec[30];

  /*
   *      Velocity statistics and parameters
   */
  real8 vAverage; /* average nodal velocity */
  real8 vStDev;   /* St.Dev of nodal velocity */

  /*
   *      I/O parameters
   */
  char dirname[MAX_STRING_LEN];
  char PrecipitatesFileName[MAX_STRING_LEN];
  char APBFileName[MAX_STRING_LEN];
  int ShearablePrecipitates;
  int PairPrecipitateShearing;
  double PrecipitateNucleationRate;
  double PrecipitateNucleationRadius;
  int PrecipitateNucleationMaxCount;

  int numIOGroups; /* number of groups into which to split tasks */
                   /* when doing parallel I/O                    */

  int savefreq, savecounter;

  real8 savepropdt, saveproptime;

  int savedensityspec[3];

  int tecplot, tecplotcounter;
  real8 tecplotdt, tecplottime;

  /*
   *      Lengths (and reciprocals) of each side of the problem
   *      space box.
   */
  real8 Lx, Ly, Lz;
  real8 invLx, invLy, invLz;

  /*
   *      General stuff
   */
  real8 springConst;
  real8 rann; /* closest distance before dislocations are */
              /* considered in contact    */

  int numBurgGroups; /* Number of groups into which different  */
                     /* burgers vectors are organized in order */
                     /* to track dislocation density by burgers*/
                     /* vector.  This number is dependent on   */
                     /* the type of mobility used.             */

  real8 disloDensity;

  real8 delSegLength; /* accumulated length of deleted segments    */
                      /* since most recent write to result_density */
                      /* property file                             */

  real8 densityChange[14]; /* For tracking density change by burgers */
                           /* vector; currently only used for BCC.   */
                           /* Values accumulated since most recent   */
                           /* write of density change data to the    */
                           /* property file.                         */

  real8 delpStrain[6], mergedelpStrain[6], totpStn[6];
  real8 delpSpin[6], mergedelpSpin[6], totpSpn[6];

  /*
   *      Added for strain decomposition and density flux decomp.
   */
  real8 totstraintensor[6];

  int imgstrgrid[6];

  /*
   *      The following two variables are used only temporarily and as such
   *      should not be specifiable in the control file, and hence, should
   *      not be bound to identifying strings via bindvar() as the other
   *      elements of this structure are.
   */
  char node_data_file[MAX_STRING_LEN];
  /* name of the file containing the    */
  /* nodal data for the run.  This data */
  /* is in the newer nodal data format, */
  /* (i.e. not contained as part of the */
  /* control file itself)               */

  /*
   *      Define the parameters used in the nodal data file
   */
  int dataFileVersion;  /* Version number of the data file */
  int nodeCount;        /* Total number of nodes in the    */
                        /* data file (all file segments)   */
  double Dimensions[3]; /* dimensions in the XYZ directions  */
  double XAxis[3];      /* local x axis coordinates  */
  double YAxis[3];      /* local y axis coordinates  */
  double ZAxis[3];      /* local z axis coordinates  */
  double GrainDiameter;
  int GeometryType;
  // 1 for cuboidal geometry, 2 for full grains and 3 for half grains
  double center[3];

  /*
   *      Define a couple factors used when calculating dislocation density.
   *      These will be calculated during initialization and will not be
   *      specified in the control parameter file.
   */
  real8 simVol;         /* Total volume of simulation */
  real8 burgVolFactor;  /* Volume factor used to convert dislocation */
                        /* length to dislocation density             */
  int Time_evo;         /* Number of stress evolution steps */
  int stress_timestep;  /* Current time step of stress evolution*/
  int GP_x, GP_y, GP_z; /* Number of grid points for stress data in a domain */
  real8 stress_Dim[3];
};

void CtrlParamInit(Param_t *param, ParamList_t *CPList);
void DataParamInit(Param_t *param, ParamList_t *DPList);
void StressParamInit(Param_t *param, ParamList_t *SPList);

#endif /* _PARAM_H */
