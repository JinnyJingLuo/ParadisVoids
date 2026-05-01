/**************************************************************************
 *
 *      Module:       Param.c
 *      Description:  Contains functions for binding control and data
 *                    file parameter names to code variables.
 *
 *      Public functions:
 *          CtrlParamInit()
 *          DataParamInit()
 *          MarkParamDisabled()
 *          MarkParamEnabled()
 *
 *************************************************************************/
#include <stdarg.h>
#include <ctype.h>
#include "Home.h"
#include "Parse.h"

/*------------------------------------------------------------------------
 *
 *      Functions:    MarkParamDisabled
 *      Description:  Explicitly mark a control file parameter
 *                    to be disabled.  Primarily used to
 *                    prevent writing to the restart file any
 *                    control file parameters that are not
 *                    applicable or appropriate to the current
 *                    execution of the code.
 *
 *      Arguments:
 *          CPList    Pointer to the parameter list structure
 *                    associated with the control file parameters.
 *          name      String containing the name of the control
 *                    parameter to be explicitly disabled.
 *
 *----------------------------------------------------------------------*/
void MarkParamDisabled(ParamList_t *CPList, char *name) {
  int paramIndex;

  if ((paramIndex = LookupParam(CPList, name)) >= 0) {
    CPList->varList[paramIndex].flags |= VFLAG_DISABLED;
  } else {
    Fatal("MarkParamDisabled: unknown parameter %s\n", name);
  }

  return;
}

/*------------------------------------------------------------------------
 *
 *      Functions:    MarkParamEnabled
 *      Description:  Explicitly mark a control file parameter
 *                    to be enabled.  Primarily used when
 *                    explicitly controlling whether a control
 *                    parameter is to be written to a restart
 *                    file.
 *
 *      Arguments:
 *          CPList    Pointer to the parameter list structure
 *                    associated with the control file parameters.
 *          name      String containing the name of the control
 *                    parameter to be explicitly enabled.
 *
 *----------------------------------------------------------------------*/
void MarkParamEnabled(ParamList_t *CPList, char *name) {
  int paramIndex;

  if ((paramIndex = LookupParam(CPList, name)) >= 0) {
    CPList->varList[paramIndex].flags &= ~VFLAG_DISABLED;
  } else {
    Fatal("MarkParamEnabled: unknown parameter %s\n", name);
  }

  return;
}

/*------------------------------------------------------------------------
 *
 *      Function:     CtrlParamInit
 *      Description:  Bind all valid control file parameters to
 *                    the associated code variables.
 *
 *      Arguments:
 *          CPList    Pointer to the parameter list structure
 *                    associated with the control file parameters.
 *
 *----------------------------------------------------------------------*/
void CtrlParamInit(Param_t *param, ParamList_t *CPList) {
  /*
   *      Note: Parameters need only be initialized if their
   *      default values are non-zero.
   */
  CPList->paramCnt = 0;
  CPList->varList = (VarData_t *)NULL;

  /*
   *      Simulation cell and processor setup
   */
  BindVar(CPList, "Simulation cell and processor setup", (void *)NULL,
          V_COMMENT, 1, VFLAG_NULL);

  BindVar(CPList, "numXdoms", &param->nXdoms, V_INT, 1, VFLAG_NULL);
  param->nXdoms = 1;

  BindVar(CPList, "numYdoms", &param->nYdoms, V_INT, 1, VFLAG_NULL);
  param->nYdoms = 1;

  BindVar(CPList, "numZdoms", &param->nZdoms, V_INT, 1, VFLAG_NULL);
  param->nZdoms = 1;

  BindVar(CPList, "numXcells", &param->nXcells, V_INT, 1, VFLAG_NULL);
  param->nXcells = 3; /* Must be >= 3 */

  BindVar(CPList, "numYcells", &param->nYcells, V_INT, 1, VFLAG_NULL);
  param->nYcells = 3; /* Must be >= 3 */

#ifdef _BGP
  BindVar(CPList, "taskMappingMode", &param->taskMappingMode, V_INT, 1,
          VFLAG_NULL);
  param->taskMappingMode = 1; /* use user-supplied decomp by default */
#endif

  BindVar(CPList, "numZcells", &param->nZcells, V_INT, 1, VFLAG_NULL);
  param->nZcells = 3; /* Must be >= 3 */

  BindVar(CPList, "BoundaryType", &param->BoundaryType, V_INT, 1, VFLAG_NULL);
  param->BoundaryType = FREE_BOUNDARY;

  BindVar(CPList, "DLBfreq", &param->DLBfreq, V_INT, 1, VFLAG_NULL);
  param->DLBfreq = 3;

  BindVar(CPList, "Ttype", &param->Ttype, V_INT, 1, VFLAG_NULL);
  param->Ttype = 0;

  BindVar(CPList, "SType", &param->SType, V_INT, 1, VFLAG_NULL);
  param->SType = 0;

  BindVar(CPList, "EnablePeierls", &param->EnablePeierls, V_INT, 1, VFLAG_NULL);
  param->EnablePeierls = 0;

  BindVar(CPList, "PeierlsScrew", &param->PeierlsScrew, V_DBL, 1, VFLAG_NULL);
  param->PeierlsScrew = 0.0;

  BindVar(CPList, "PeierlsEdge", &param->PeierlsEdge, V_DBL, 1, VFLAG_NULL);
  param->PeierlsEdge = 0.0;

  BindVar(CPList, "PeierlsOthers", &param->PeierlsOthers, V_DBL, 1, VFLAG_NULL);
  param->PeierlsOthers = 0.0;

  BindVar(CPList, "EnableTwinPlaneCrossSlip", &param->EnableTwinPlaneCrossSlip,
          V_INT, 1, VFLAG_NULL);
  param->EnableTwinPlaneCrossSlip = 0;
  BindVar(CPList, "A", &param->A, V_DBL, 1, VFLAG_NULL);
  param->A = 1.0;
  BindVar(CPList, "B", &param->B, V_DBL, 1, VFLAG_NULL);
  param->B = 0.0;
  BindVar(CPList, "C", &param->C, V_DBL, 1, VFLAG_NULL);
  param->C = 0.0;
  BindVar(CPList, "X_1", &param->X_1, V_DBL, 1, VFLAG_NULL);
  param->X_1 = 0.0;
  BindVar(CPList, "Y_1", &param->Y_1, V_DBL, 1, VFLAG_NULL);
  param->Y_1 = 0.0;
  BindVar(CPList, "Z_1", &param->Z_1, V_DBL, 1, VFLAG_NULL);
  param->Z_1 = 0.0;

  /*
   *      Simulation time and timestepping controls
   */
  BindVar(CPList, "Simulation time and timestepping controls", (void *)NULL,
          V_COMMENT, 1, VFLAG_NULL);

  BindVar(CPList, "cycleStart", &param->cycleStart, V_INT, 1, VFLAG_NULL);
  param->cycleStart = 1;

  BindVar(CPList, "maxstep", &param->maxstep, V_INT, 1, VFLAG_NULL);
  param->maxstep = 100;

  BindVar(CPList, "timeNow", &param->timeNow, V_DBL, 1, VFLAG_NULL);
  param->timeNow = 0.0;

  BindVar(CPList, "deltaTT", &param->deltaTT, V_DBL, 1, VFLAG_NULL);
  param->deltaTT = 1.0E-11;

  BindVar(CPList, "rTol", &param->rTol, V_DBL, 1, VFLAG_NULL);
  param->rTol = -1.0;

  BindVar(CPList, "rmax", &param->rmax, V_DBL, 1, VFLAG_NULL);
  param->rmax = 100;

  /*
   *      Discretization controls and controls for topological changes
   */
  BindVar(CPList, "Discretization and topological change controls",
          (void *)NULL, V_COMMENT, 1, VFLAG_NULL);

  BindVar(CPList, "maxSeg", &param->maxSeg, V_DBL, 1, VFLAG_NULL);
  param->maxSeg = -1.0;

  BindVar(CPList, "minSeg", &param->minSeg, V_DBL, 1, VFLAG_NULL);
  param->minSeg = -1.0;

  BindVar(CPList, "splitMultiNodeFreq", &param->splitMultiNodeFreq, V_INT, 1,
          VFLAG_NULL);
  param->splitMultiNodeFreq = 1;

  BindVar(CPList, "collisionMethod", &param->collisionMethod, V_INT, 1,
          VFLAG_NULL);
  param->collisionMethod = 2;

  /*
   *      Identify tables needed for remote force calculations if
   *      FMM is not enabled
   */
  BindVar(CPList, "Tables for non-FMM far-field force calcs", (void *)NULL,
          V_COMMENT, 1, VFLAG_NULL);

  BindVar(CPList, "Rijmfile", param->Rijmfile, V_STRING, 1, VFLAG_NULL);
  strcpy(param->Rijmfile, "inputs/Rijm.cube.out");

  BindVar(CPList, "RijmPBCfile", param->RijmPBCfile, V_STRING, 1, VFLAG_NULL);
  strcpy(param->RijmPBCfile, "inputs/RijmPBC.cube.out");

  /*
   *      Loading condition parameters
   */
  BindVar(CPList, "Loading conditions", (void *)NULL, V_COMMENT, 1, VFLAG_NULL);

  BindVar(CPList, "TempK", &param->TempK, V_DBL, 1, VFLAG_NULL);
  param->TempK = 600;

  BindVar(CPList, "loadType", &param->loadType, V_INT, 1, VFLAG_NULL);

  BindVar(CPList, "appliedStress", param->appliedStress, V_DBL, 6, VFLAG_NULL);

  BindVar(CPList, "eRate", &param->eRate, V_DBL, 1, VFLAG_NULL);
  param->eRate = 1.0;

  BindVar(CPList, "indxErate", &param->indxErate, V_INT, 1, VFLAG_NULL);
  param->indxErate = 1;

  BindVar(CPList, "edotdir", param->edotdir, V_DBL, 3, VFLAG_NULL);
  param->edotdir[0] = 1.0;
  param->edotdir[1] = 0.0;
  param->edotdir[2] = 0.0;

  BindVar(CPList, "cTimeOld", &param->cTimeOld, V_DBL, 1, VFLAG_NULL);

  BindVar(CPList, "dCyclicStrain", &param->dCyclicStrain, V_DBL, 1, VFLAG_NULL);

  BindVar(CPList, "netCyclicStrain", &param->netCyclicStrain, V_DBL, 1,
          VFLAG_NULL);

  BindVar(CPList, "numLoadCycle", &param->numLoadCycle, V_INT, 1, VFLAG_NULL);

  BindVar(CPList, "eAmp", &param->eAmp, V_DBL, 1, VFLAG_NULL);

  // Added by Ahmed M. Hussein, July 30th 2012
  BindVar(CPList, "HandleCrossSlip", &param->HandleCrossSlip, V_INT, 1,
          VFLAG_NULL);
  param->HandleCrossSlip = 1;
  /*BindVar(CPList, "Handle",&param->Handle,V_INT,1,VFLAG_NULL);
  param->Handle = 0;*/
  BindVar(CPList, "CrossSlipHandlingFrequency",
          &param->CrossSlipHandlingFrequency, V_INT, 1, VFLAG_NULL);
  param->CrossSlipHandlingFrequency = 1000;

  BindVar(CPList, "BulkCrossSlipActivationEnergy",
          &param->BulkCrossSlipActivationEnergy, V_DBL, 1, VFLAG_NULL);
  param->BulkCrossSlipActivationEnergy = 0.8 * 1.60217646E-19;
  BindVar(CPList, "SurfaceCrossSlipActivationEnergy",
          &param->SurfaceCrossSlipActivationEnergy, V_DBL, 1, VFLAG_NULL);
  param->SurfaceCrossSlipActivationEnergy = 0.2 * 1.60217646E-19;
  BindVar(CPList, "HirthCrossSlipActivationEnergy",
          &param->HirthCrossSlipActivationEnergy, V_DBL, 1, VFLAG_NULL);
  param->HirthCrossSlipActivationEnergy = 0.2 * 1.60217646E-19;
  BindVar(CPList, "LCLockCrossSlipActivationEnergy",
          &param->LCLockCrossSlipActivationEnergy, V_DBL, 1, VFLAG_NULL);
  param->LCLockCrossSlipActivationEnergy = 0.6 * 1.60217646E-19;
  BindVar(CPList, "GlideLockCrossSlipActivationEnergy",
          &param->GlideLockCrossSlipActivationEnergy, V_DBL, 1, VFLAG_NULL);
  param->GlideLockCrossSlipActivationEnergy = 0.5 * 1.60217646E-19;

  BindVar(CPList, "BulkCrossSlipActivationVolumeFactor",
          &param->BulkCrossSlipActivationVolumeFactor, V_DBL, 1, VFLAG_NULL);
  param->BulkCrossSlipActivationVolumeFactor = 20.0;
  BindVar(CPList, "SurfaceCrossSlipActivationVolumeFactor",
          &param->SurfaceCrossSlipActivationVolumeFactor, V_DBL, 1, VFLAG_NULL);
  param->SurfaceCrossSlipActivationVolumeFactor = 20.0;
  BindVar(CPList, "HirthCrossSlipActivationVolumeFactor",
          &param->HirthCrossSlipActivationVolumeFactor, V_DBL, 1, VFLAG_NULL);
  param->HirthCrossSlipActivationVolumeFactor = 20.0;
  BindVar(CPList, "LCLockCrossSlipActivationVolumeFactor",
          &param->LCLockCrossSlipActivationVolumeFactor, V_DBL, 1, VFLAG_NULL);
  param->LCLockCrossSlipActivationVolumeFactor = 20.0;
  BindVar(CPList, "GlideLockCrossSlipActivationVolumeFactor",
          &param->GlideLockCrossSlipActivationVolumeFactor, V_DBL, 1,
          VFLAG_NULL);
  param->GlideLockCrossSlipActivationVolumeFactor = 20.0;

  BindVar(CPList, "BulkCrossSlipReferenceLength",
          &param->BulkCrossSlipReferenceLength, V_DBL, 1, VFLAG_NULL);
  param->BulkCrossSlipReferenceLength = 4.0E3;
  BindVar(CPList, "SurfaceCrossSlipReferenceLength",
          &param->SurfaceCrossSlipReferenceLength, V_DBL, 1, VFLAG_NULL);
  param->SurfaceCrossSlipReferenceLength = 4.0E3;
  BindVar(CPList, "HirthCrossSlipReferenceLength",
          &param->HirthCrossSlipReferenceLength, V_DBL, 1, VFLAG_NULL);
  param->HirthCrossSlipReferenceLength = 4.0E3;
  BindVar(CPList, "LCLockCrossSlipReferenceLength",
          &param->LCLockCrossSlipReferenceLength, V_DBL, 1, VFLAG_NULL);
  param->LCLockCrossSlipReferenceLength = 4.0E3;
  BindVar(CPList, "GlideLockCrossSlipReferenceLength",
          &param->GlideLockCrossSlipReferenceLength, V_DBL, 1, VFLAG_NULL);
  param->GlideLockCrossSlipReferenceLength = 4.0E3;

  BindVar(CPList, "SurfaceCrossSlipLength", &param->SurfaceCrossSlipLength,
          V_DBL, 1, VFLAG_NULL);
  param->SurfaceCrossSlipLength = 5.0;
  BindVar(CPList, "HirthCrossSlipLength", &param->HirthCrossSlipLength, V_DBL,
          1, VFLAG_NULL);
  param->HirthCrossSlipLength = 10.0;
  BindVar(CPList, "LCLockCrossSlipLength", &param->LCLockCrossSlipLength, V_DBL,
          1, VFLAG_NULL);
  param->LCLockCrossSlipLength = 10.0;
  BindVar(CPList, "GlideLockCrossSlipLength", &param->GlideLockCrossSlipLength,
          V_DBL, 1, VFLAG_NULL);
  param->GlideLockCrossSlipLength = 10.0;

  BindVar(CPList, "BulkCrossSlipFrequency", &param->BulkCrossSlipFrequency,
          V_DBL, 1, VFLAG_NULL);
  param->BulkCrossSlipFrequency = 5.0E17;
  BindVar(CPList, "SurfaceCrossSlipFrequency",
          &param->SurfaceCrossSlipFrequency, V_DBL, 1, VFLAG_NULL);
  param->SurfaceCrossSlipFrequency = 5.0E17;
  BindVar(CPList, "HirthCrossSlipFrequency", &param->HirthCrossSlipFrequency,
          V_DBL, 1, VFLAG_NULL);
  param->HirthCrossSlipFrequency = 5.0E17;
  BindVar(CPList, "LCLockCrossSlipFrequency", &param->LCLockCrossSlipFrequency,
          V_DBL, 1, VFLAG_NULL);
  param->LCLockCrossSlipFrequency = 5.0E17;
  BindVar(CPList, "GlideLockCrossSlipFrequency",
          &param->GlideLockCrossSlipFrequency, V_DBL, 1, VFLAG_NULL);
  param->GlideLockCrossSlipFrequency = 5.0E17;

  BindVar(CPList, "BulkStackingFaultEnergy", &param->BulkStackingFaultEnergy,
          V_DBL, 1, VFLAG_NULL);
  param->BulkStackingFaultEnergy = 0.120;
  BindVar(CPList, "SurfaceStackingFaultEnergy",
          &param->SurfaceStackingFaultEnergy, V_DBL, 1, VFLAG_NULL);
  param->SurfaceStackingFaultEnergy = 0.120;
  BindVar(CPList, "HirthStackingFaultEnergy", &param->HirthStackingFaultEnergy,
          V_DBL, 1, VFLAG_NULL);
  param->HirthStackingFaultEnergy = 0.120;
  BindVar(CPList, "LCLockStackingFaultEnergy",
          &param->LCLockStackingFaultEnergy, V_DBL, 1, VFLAG_NULL);
  param->LCLockStackingFaultEnergy = 0.120;
  BindVar(CPList, "GlideLockStackingFaultEnergy",
          &param->GlideLockStackingFaultEnergy, V_DBL, 1, VFLAG_NULL);
  param->GlideLockStackingFaultEnergy = 0.120;

  BindVar(CPList, "BulkDoverBRatio", &param->BulkDoverBRatio, V_DBL, 1,
          VFLAG_NULL);
  param->BulkDoverBRatio = 5.0;
  BindVar(CPList, "SurfaceDoverBRatio", &param->SurfaceDoverBRatio, V_DBL, 1,
          VFLAG_NULL);
  param->SurfaceDoverBRatio = 5.0;
  BindVar(CPList, "HirthDoverBRatio", &param->HirthDoverBRatio, V_DBL, 1,
          VFLAG_NULL);
  param->HirthDoverBRatio = 5.0;
  BindVar(CPList, "LCLockDoverBRatio", &param->LCLockDoverBRatio, V_DBL, 1,
          VFLAG_NULL);
  param->LCLockDoverBRatio = 5.0;
  BindVar(CPList, "GlideLockDoverBRatio", &param->GlideLockDoverBRatio, V_DBL,
          1, VFLAG_NULL);
  param->GlideLockDoverBRatio = 5.0;

  // cross slip statistics
  BindVar(CPList, "BulkCrossSlipEventsCount", &param->BulkCrossSlipEventsCount,
          V_INT, 1, VFLAG_NULL);
  param->BulkCrossSlipEventsCount = 0;
  BindVar(CPList, "SurfaceCrossSlipEventsCount",
          &param->SurfaceCrossSlipEventsCount, V_INT, 1, VFLAG_NULL);
  param->SurfaceCrossSlipEventsCount = 0;
  BindVar(CPList, "AttractiveCrossSlipEventsCount",
          &param->AttractiveCrossSlipEventsCount, V_INT, 1, VFLAG_NULL);
  param->AttractiveCrossSlipEventsCount = 0;
  BindVar(CPList, "RepulsiveCrossSlipEventsCount",
          &param->RepulsiveCrossSlipEventsCount, V_INT, 1, VFLAG_NULL);
  param->RepulsiveCrossSlipEventsCount = 0;

  BindVar(CPList, "TotalBulkCrossSlippedChainsLength",
          &param->TotalBulkCrossSlippedChainsLength, V_DBL, 1, VFLAG_NULL);
  param->TotalBulkCrossSlippedChainsLength = 0.0;
  BindVar(CPList, "TotalSurfaceCrossSlippedChainsLength",
          &param->TotalSurfaceCrossSlippedChainsLength, V_DBL, 1, VFLAG_NULL);
  param->TotalSurfaceCrossSlippedChainsLength = 0.0;
  BindVar(CPList, "TotalAttractiveCrossSlippedChainsLength",
          &param->TotalAttractiveCrossSlippedChainsLength, V_DBL, 1,
          VFLAG_NULL);
  param->TotalAttractiveCrossSlippedChainsLength = 0.0;
  BindVar(CPList, "TotalRepulsiveCrossSlippedChainsLength",
          &param->TotalRepulsiveCrossSlippedChainsLength, V_DBL, 1, VFLAG_NULL);
  param->TotalRepulsiveCrossSlippedChainsLength = 0.0;

  // load parameters for FEM model
  BindVar(CPList, "CyclicLoadAmplitude", &param->CyclicLoadAmplitude, V_DBL, 1,
          VFLAG_NULL);
  param->CyclicLoadAmplitude = 0.0;
  BindVar(CPList, "CyclicLoadFrequency", &param->CyclicLoadFrequency, V_DBL, 1,
          VFLAG_NULL);
  param->CyclicLoadFrequency = 0.0;
  BindVar(CPList, "InitialTemperature", &param->InitialTemperature, V_DBL, 1,
          VFLAG_NULL);
  param->InitialTemperature = 0.0;
  BindVar(CPList, "HeatingRate", &param->HeatingRate, V_DBL, 1, VFLAG_NULL);
  param->HeatingRate = 0.0;

  // material properties for FEM model, default is Nickel
  BindVar(CPList, "YoungsModulus", &param->YoungsModulus, V_DBL, 1, VFLAG_NULL);
  param->YoungsModulus = 2.0E11;
  BindVar(CPList, "PoissonsRatio", &param->PoissonsRatio, V_DBL, 1, VFLAG_NULL);
  param->PoissonsRatio = 0.31;
  BindVar(CPList, "MassDensity", &param->MassDensity, V_DBL, 1, VFLAG_NULL);
  param->MassDensity = 8912.0;
  BindVar(CPList, "ThermalConductivity", &param->ThermalConductivity, V_DBL, 1,
          VFLAG_NULL);
  param->ThermalConductivity = 90.9;
  BindVar(CPList, "ThermalExpansionCoefficient",
          &param->ThermalExpansionCoefficient, V_DBL, 1, VFLAG_NULL);
  param->ThermalExpansionCoefficient = 1.3E-6;
  BindVar(CPList, "SpecificHeatCapacity", &param->SpecificHeatCapacity, V_DBL,
          1, VFLAG_NULL);
  param->SpecificHeatCapacity = 440.0;
  BindVar(CPList, "ReferenceTemperature", &param->ReferenceTemperature, V_DBL,
          1, VFLAG_NULL);
  param->ReferenceTemperature = 298.0;
  BindVar(CPList, "FEMTimeStep", &param->FEMTimeStep, V_DBL, 1, VFLAG_NULL);
  param->FEMTimeStep = 1.0E-7;

  BindVar(CPList, "APBEnrgy", &param->APBEnrgy, V_DBL, 1, VFLAG_NULL);
  param->APBEnrgy = 0.1;

  BindVar(CPList, "FrictionStress", &param->FrictionStress, V_DBL, 1,
          VFLAG_NULL);
  param->FrictionStress = -1.0;

  BindVar(CPList, "SolidSolutionStrength", &param->SolidSolutionStrength, V_DBL,
          1, VFLAG_NULL);
  param->SolidSolutionStrength = -1.0;

  BindVar(CPList, "KWStress", &param->KWStress, V_DBL, 1, VFLAG_NULL);
  param->KWStress = -1.0;

  BindVar(CPList, "voidR", &param->voidR, V_DBL, 1, VFLAG_NULL);
  param->voidR = -1.0;

  BindVar(CPList, "VelocityDampingSteps", &param->VelocityDampingSteps, V_INT,
          1, VFLAG_NULL);
  param->VelocityDampingSteps = 6;

  BindVar(CPList, "EnableFEM", &param->EnableFEM, V_INT, 1, VFLAG_NULL);
  param->EnableFEM = 0;

  BindVar(CPList, "FEMCorrectionFrequency", &param->FEMCorrectionFrequency,
          V_INT, 1, VFLAG_NULL);
  param->FEMCorrectionFrequency = 10000;

  BindVar(CPList, "LoopExtractionFrequency", &param->LoopExtractionFrequency,
          V_INT, 1, VFLAG_NULL);
  param->LoopExtractionFrequency = 0;

  BindVar(CPList, "FEMRestartFileName", param->FEMRestartFileName, V_STRING, 1,
          VFLAG_NULL);
  strcpy(param->FEMRestartFileName, "");

  BindVar(CPList, "AnnihilateSegments", &param->AnnihilateSegments, V_INT, 1,
          VFLAG_NULL);
  param->AnnihilateSegments = 0;

  BindVar(CPList, "AnnihilationDistance", &param->AnnihilationDistance, V_DBL,
          1, VFLAG_NULL);
  param->AnnihilationDistance = 6.0;

  /*
   *      To specify input sample axes in laboratory frame
   */
  BindVar(CPList, "useLabFrame", &param->useLabFrame, V_INT, 1, VFLAG_NULL);

  BindVar(CPList, "labFrameXDir", param->labFrameXDir, V_DBL, 3, VFLAG_NULL);
  param->labFrameXDir[0] = 1.0;
  param->labFrameXDir[1] = 0.0;
  param->labFrameXDir[2] = 0.0;

  BindVar(CPList, "labFrameYDir", param->labFrameYDir, V_DBL, 3, VFLAG_NULL);
  param->labFrameYDir[0] = 0.0;
  param->labFrameYDir[1] = 1.0;
  param->labFrameYDir[2] = 0.0;

  BindVar(CPList, "labFrameZDir", param->labFrameZDir, V_DBL, 3, VFLAG_NULL);
  param->labFrameZDir[0] = 1.0;
  param->labFrameZDir[1] = 0.0;
  param->labFrameZDir[2] = 0.0;

  /*
   *	Cross-slip per slip system
   */
  BindVar(CPList, "Cross-slip per slip system", (void *)NULL, V_COMMENT, 1,
          VFLAG_NULL);
  /*BindVar(CPList,
  "HandleCrossSlipPerSlipSystem",&param->HandleCrossSlipPerSlipSystem,V_INT,1,VFLAG_NULL);
  param->HandleCrossSlipPerSlipSystem = 1;*/
  BindVar(CPList, "HandleCrossSlipPerSlipSystem",
          &param->HandleCrossSlipPerSlipSystem, V_INT, 1, VFLAG_NULL);
  param->HandleCrossSlipPerSlipSystem = 1;
  BindVar(CPList, "SlipSystem1BulkCrossSlipEventsCount",
          &param->SlipSystem1BulkCrossSlipEventsCount, V_INT, 1, VFLAG_NULL);
  param->SlipSystem1BulkCrossSlipEventsCount = 0;
  BindVar(CPList, "SlipSystem2BulkCrossSlipEventsCount",
          &param->SlipSystem2BulkCrossSlipEventsCount, V_INT, 1, VFLAG_NULL);
  param->SlipSystem2BulkCrossSlipEventsCount = 0;
  BindVar(CPList, "SlipSystem3BulkCrossSlipEventsCount",
          &param->SlipSystem3BulkCrossSlipEventsCount, V_INT, 1, VFLAG_NULL);
  param->SlipSystem3BulkCrossSlipEventsCount = 0;
  BindVar(CPList, "SlipSystem4BulkCrossSlipEventsCount",
          &param->SlipSystem4BulkCrossSlipEventsCount, V_INT, 1, VFLAG_NULL);
  param->SlipSystem4BulkCrossSlipEventsCount = 0;
  BindVar(CPList, "SlipSystem5BulkCrossSlipEventsCount",
          &param->SlipSystem5BulkCrossSlipEventsCount, V_INT, 1, VFLAG_NULL);
  param->SlipSystem5BulkCrossSlipEventsCount = 0;
  BindVar(CPList, "SlipSystem6BulkCrossSlipEventsCount",
          &param->SlipSystem6BulkCrossSlipEventsCount, V_INT, 1, VFLAG_NULL);
  param->SlipSystem6BulkCrossSlipEventsCount = 0;
  BindVar(CPList, "SlipSystem7BulkCrossSlipEventsCount",
          &param->SlipSystem7BulkCrossSlipEventsCount, V_INT, 1, VFLAG_NULL);
  param->SlipSystem7BulkCrossSlipEventsCount = 0;
  BindVar(CPList, "SlipSystem8BulkCrossSlipEventsCount",
          &param->SlipSystem8BulkCrossSlipEventsCount, V_INT, 1, VFLAG_NULL);
  param->SlipSystem8BulkCrossSlipEventsCount = 0;
  BindVar(CPList, "SlipSystem9BulkCrossSlipEventsCount",
          &param->SlipSystem9BulkCrossSlipEventsCount, V_INT, 1, VFLAG_NULL);
  param->SlipSystem9BulkCrossSlipEventsCount = 0;
  BindVar(CPList, "SlipSystem10BulkCrossSlipEventsCount",
          &param->SlipSystem10BulkCrossSlipEventsCount, V_INT, 1, VFLAG_NULL);
  param->SlipSystem10BulkCrossSlipEventsCount = 0;
  BindVar(CPList, "SlipSystem11BulkCrossSlipEventsCount",
          &param->SlipSystem11BulkCrossSlipEventsCount, V_INT, 1, VFLAG_NULL);
  param->SlipSystem11BulkCrossSlipEventsCount = 0;
  BindVar(CPList, "SlipSystem12BulkCrossSlipEventsCount",
          &param->SlipSystem12BulkCrossSlipEventsCount, V_INT, 1, VFLAG_NULL);
  param->SlipSystem12BulkCrossSlipEventsCount = 0;
  BindVar(CPList, "OtherSlipSystemBulkCrossSlipEventsCount",
          &param->OtherSlipSystemBulkCrossSlipEventsCount, V_INT, 1,
          VFLAG_NULL);
  param->OtherSlipSystemBulkCrossSlipEventsCount = 0;
  BindVar(CPList, "SlipSystem1SurfaceCrossSlipEventsCount",
          &param->SlipSystem1SurfaceCrossSlipEventsCount, V_INT, 1, VFLAG_NULL);
  param->SlipSystem1SurfaceCrossSlipEventsCount = 0;
  BindVar(CPList, "SlipSystem2SurfaceCrossSlipEventsCount",
          &param->SlipSystem2SurfaceCrossSlipEventsCount, V_INT, 1, VFLAG_NULL);
  param->SlipSystem2SurfaceCrossSlipEventsCount = 0;
  BindVar(CPList, "SlipSystem3SurfaceCrossSlipEventsCount",
          &param->SlipSystem3SurfaceCrossSlipEventsCount, V_INT, 1, VFLAG_NULL);
  param->SlipSystem3SurfaceCrossSlipEventsCount = 0;
  BindVar(CPList, "SlipSystem4SurfaceCrossSlipEventsCount",
          &param->SlipSystem4SurfaceCrossSlipEventsCount, V_INT, 1, VFLAG_NULL);
  param->SlipSystem4SurfaceCrossSlipEventsCount = 0;
  BindVar(CPList, "SlipSystem5SurfaceCrossSlipEventsCount",
          &param->SlipSystem5SurfaceCrossSlipEventsCount, V_INT, 1, VFLAG_NULL);
  param->SlipSystem5SurfaceCrossSlipEventsCount = 0;
  BindVar(CPList, "SlipSystem6SurfaceCrossSlipEventsCount",
          &param->SlipSystem6SurfaceCrossSlipEventsCount, V_INT, 1, VFLAG_NULL);
  param->SlipSystem6SurfaceCrossSlipEventsCount = 0;
  BindVar(CPList, "SlipSystem7SurfaceCrossSlipEventsCount",
          &param->SlipSystem7SurfaceCrossSlipEventsCount, V_INT, 1, VFLAG_NULL);
  param->SlipSystem7SurfaceCrossSlipEventsCount = 0;
  BindVar(CPList, "SlipSystem8SurfaceCrossSlipEventsCount",
          &param->SlipSystem8SurfaceCrossSlipEventsCount, V_INT, 1, VFLAG_NULL);
  param->SlipSystem8SurfaceCrossSlipEventsCount = 0;
  BindVar(CPList, "SlipSystem9SurfaceCrossSlipEventsCount",
          &param->SlipSystem9SurfaceCrossSlipEventsCount, V_INT, 1, VFLAG_NULL);
  param->SlipSystem9SurfaceCrossSlipEventsCount = 0;
  BindVar(CPList, "SlipSystem10SurfaceCrossSlipEventsCount",
          &param->SlipSystem10SurfaceCrossSlipEventsCount, V_INT, 1,
          VFLAG_NULL);
  param->SlipSystem10SurfaceCrossSlipEventsCount = 0;
  BindVar(CPList, "SlipSystem11SurfaceCrossSlipEventsCount",
          &param->SlipSystem11SurfaceCrossSlipEventsCount, V_INT, 1,
          VFLAG_NULL);
  param->SlipSystem11SurfaceCrossSlipEventsCount = 0;
  BindVar(CPList, "SlipSystem12SurfaceCrossSlipEventsCount",
          &param->SlipSystem12SurfaceCrossSlipEventsCount, V_INT, 1,
          VFLAG_NULL);
  param->SlipSystem12SurfaceCrossSlipEventsCount = 0;
  BindVar(CPList, "OtherSlipSystemSurfaceCrossSlipEventsCount",
          &param->OtherSlipSystemSurfaceCrossSlipEventsCount, V_INT, 1,
          VFLAG_NULL);
  param->OtherSlipSystemSurfaceCrossSlipEventsCount = 0;
  BindVar(CPList, "SlipSystem1AttractiveCrossSlipEventsCount",
          &param->SlipSystem1AttractiveCrossSlipEventsCount, V_INT, 1,
          VFLAG_NULL);
  param->SlipSystem1AttractiveCrossSlipEventsCount = 0;
  BindVar(CPList, "SlipSystem2AttractiveCrossSlipEventsCount",
          &param->SlipSystem2AttractiveCrossSlipEventsCount, V_INT, 1,
          VFLAG_NULL);
  param->SlipSystem2AttractiveCrossSlipEventsCount = 0;
  BindVar(CPList, "SlipSystem3AttractiveCrossSlipEventsCount",
          &param->SlipSystem3AttractiveCrossSlipEventsCount, V_INT, 1,
          VFLAG_NULL);
  param->SlipSystem3AttractiveCrossSlipEventsCount = 0;
  BindVar(CPList, "SlipSystem4AttractiveCrossSlipEventsCount",
          &param->SlipSystem4AttractiveCrossSlipEventsCount, V_INT, 1,
          VFLAG_NULL);
  param->SlipSystem4AttractiveCrossSlipEventsCount = 0;
  BindVar(CPList, "SlipSystem5AttractiveCrossSlipEventsCount",
          &param->SlipSystem5AttractiveCrossSlipEventsCount, V_INT, 1,
          VFLAG_NULL);
  param->SlipSystem5AttractiveCrossSlipEventsCount = 0;
  BindVar(CPList, "SlipSystem6AttractiveCrossSlipEventsCount",
          &param->SlipSystem6AttractiveCrossSlipEventsCount, V_INT, 1,
          VFLAG_NULL);
  param->SlipSystem6AttractiveCrossSlipEventsCount = 0;
  BindVar(CPList, "SlipSystem7AttractiveCrossSlipEventsCount",
          &param->SlipSystem7AttractiveCrossSlipEventsCount, V_INT, 1,
          VFLAG_NULL);
  param->SlipSystem7AttractiveCrossSlipEventsCount = 0;
  BindVar(CPList, "SlipSystem8AttractiveCrossSlipEventsCount",
          &param->SlipSystem8AttractiveCrossSlipEventsCount, V_INT, 1,
          VFLAG_NULL);
  param->SlipSystem8AttractiveCrossSlipEventsCount = 0;
  BindVar(CPList, "SlipSystem9AttractiveCrossSlipEventsCount",
          &param->SlipSystem9AttractiveCrossSlipEventsCount, V_INT, 1,
          VFLAG_NULL);
  param->SlipSystem9AttractiveCrossSlipEventsCount = 0;
  BindVar(CPList, "SlipSystem10AttractiveCrossSlipEventsCount",
          &param->SlipSystem10AttractiveCrossSlipEventsCount, V_INT, 1,
          VFLAG_NULL);
  param->SlipSystem10AttractiveCrossSlipEventsCount = 0;
  BindVar(CPList, "SlipSystem11AttractiveCrossSlipEventsCount",
          &param->SlipSystem11AttractiveCrossSlipEventsCount, V_INT, 1,
          VFLAG_NULL);
  param->SlipSystem11AttractiveCrossSlipEventsCount = 0;
  BindVar(CPList, "SlipSystem12AttractiveCrossSlipEventsCount",
          &param->SlipSystem12AttractiveCrossSlipEventsCount, V_INT, 1,
          VFLAG_NULL);
  param->SlipSystem12AttractiveCrossSlipEventsCount = 0;
  BindVar(CPList, "OtherSlipSystemAttractiveCrossSlipEventsCount",
          &param->OtherSlipSystemAttractiveCrossSlipEventsCount, V_INT, 1,
          VFLAG_NULL);
  param->OtherSlipSystemAttractiveCrossSlipEventsCount = 0;
  BindVar(CPList, "SlipSystem1RepulsiveCrossSlipEventsCount",
          &param->SlipSystem1RepulsiveCrossSlipEventsCount, V_INT, 1,
          VFLAG_NULL);
  param->SlipSystem1RepulsiveCrossSlipEventsCount = 0;
  BindVar(CPList, "SlipSystem2RepulsiveCrossSlipEventsCount",
          &param->SlipSystem2RepulsiveCrossSlipEventsCount, V_INT, 1,
          VFLAG_NULL);
  param->SlipSystem2RepulsiveCrossSlipEventsCount = 0;
  BindVar(CPList, "SlipSystem3RepulsiveCrossSlipEventsCount",
          &param->SlipSystem3RepulsiveCrossSlipEventsCount, V_INT, 1,
          VFLAG_NULL);
  param->SlipSystem3RepulsiveCrossSlipEventsCount = 0;
  BindVar(CPList, "SlipSystem4RepulsiveCrossSlipEventsCount",
          &param->SlipSystem4RepulsiveCrossSlipEventsCount, V_INT, 1,
          VFLAG_NULL);
  param->SlipSystem4RepulsiveCrossSlipEventsCount = 0;
  BindVar(CPList, "SlipSystem5RepulsiveCrossSlipEventsCount",
          &param->SlipSystem5RepulsiveCrossSlipEventsCount, V_INT, 1,
          VFLAG_NULL);
  param->SlipSystem5RepulsiveCrossSlipEventsCount = 0;
  BindVar(CPList, "SlipSystem6RepulsiveCrossSlipEventsCount",
          &param->SlipSystem6RepulsiveCrossSlipEventsCount, V_INT, 1,
          VFLAG_NULL);
  param->SlipSystem6RepulsiveCrossSlipEventsCount = 0;
  BindVar(CPList, "SlipSystem7RepulsiveCrossSlipEventsCount",
          &param->SlipSystem7RepulsiveCrossSlipEventsCount, V_INT, 1,
          VFLAG_NULL);
  param->SlipSystem7RepulsiveCrossSlipEventsCount = 0;
  BindVar(CPList, "SlipSystem8RepulsiveCrossSlipEventsCount",
          &param->SlipSystem8RepulsiveCrossSlipEventsCount, V_INT, 1,
          VFLAG_NULL);
  param->SlipSystem8RepulsiveCrossSlipEventsCount = 0;
  BindVar(CPList, "SlipSystem9RepulsiveCrossSlipEventsCount",
          &param->SlipSystem9RepulsiveCrossSlipEventsCount, V_INT, 1,
          VFLAG_NULL);
  param->SlipSystem9RepulsiveCrossSlipEventsCount = 0;
  BindVar(CPList, "SlipSystem10RepulsiveCrossSlipEventsCount",
          &param->SlipSystem10RepulsiveCrossSlipEventsCount, V_INT, 1,
          VFLAG_NULL);
  param->SlipSystem10RepulsiveCrossSlipEventsCount = 0;
  BindVar(CPList, "SlipSystem11RepulsiveCrossSlipEventsCount",
          &param->SlipSystem11RepulsiveCrossSlipEventsCount, V_INT, 1,
          VFLAG_NULL);
  param->SlipSystem11RepulsiveCrossSlipEventsCount = 0;
  BindVar(CPList, "SlipSystem12RepulsiveCrossSlipEventsCount",
          &param->SlipSystem12RepulsiveCrossSlipEventsCount, V_INT, 1,
          VFLAG_NULL);
  param->SlipSystem12RepulsiveCrossSlipEventsCount = 0;
  BindVar(CPList, "OtherSlipSystemRepulsiveCrossSlipEventsCount",
          &param->OtherSlipSystemRepulsiveCrossSlipEventsCount, V_INT, 1,
          VFLAG_NULL);
  param->OtherSlipSystemRepulsiveCrossSlipEventsCount = 0;

  /*
   *      Parameters for material specific constants and mobility values
   */
  BindVar(CPList, "Material and mobility parameters", (void *)NULL, V_COMMENT,
          1, VFLAG_NULL);

  BindVar(CPList, "mobilityLaw", param->mobilityLaw, V_STRING, 1, VFLAG_NULL);
  strcpy(param->mobilityLaw, "BCC_0");

  BindVar(CPList, "vacancyConc", &param->vacancyConc, V_DBL, 1, VFLAG_NULL);
  param->vacancyConc = -1.0; /* signifies not used */

  BindVar(CPList, "vacancyConcEquilibrium", &param->vacancyConcEquilibrium,
          V_DBL, 1, VFLAG_NULL);
  param->vacancyConcEquilibrium = -1.0; /* signifies not used */

  BindVar(CPList, "shearModulus", &param->shearModulus, V_DBL, 1, VFLAG_NULL);
  param->shearModulus = 6.488424e+10; /* Ta: units=Pa @t=600K, p=0Pa*/

  BindVar(CPList, "pois", &param->pois, V_DBL, 1, VFLAG_NULL);
  param->pois = 3.327533e-01; /* Ta: temp=600K, pressure=0Pa */

  BindVar(CPList, "burgMag", &param->burgMag, V_DBL, 1, VFLAG_NULL);
  param->burgMag = 2.875401e-10; /* Ta: units=m @t=600K, pressure=0Pa */

  BindVar(CPList, "rc", &param->rc, V_DBL, 1, VFLAG_NULL);
  param->rc = -1.0;

  BindVar(CPList, "Ecore", &param->Ecore, V_DBL, 1, VFLAG_NULL);
  param->Ecore = -1.0;

  BindVar(CPList, "MobScrew", &param->MobScrew, V_DBL, 1, VFLAG_NULL);
  param->MobScrew = 10.0;

  BindVar(CPList, "MobEdge", &param->MobEdge, V_DBL, 1, VFLAG_NULL);
  param->MobEdge = 10.0;

  BindVar(CPList, "MobClimb", &param->MobClimb, V_DBL, 1, VFLAG_NULL);
  param->MobClimb = 1.0e-02;

  BindVar(CPList, "MobM", &param->MobM, V_DBL, 1, VFLAG_NULL);
  param->MobM = 0.1;

  BindVar(CPList, "MobB0", &param->MobB0, V_DBL, 1, VFLAG_NULL);
  param->MobB0 = 0.1;
  /*
   *      List of burgers vectors/line directions to be considered sessile
   */
  BindVar(CPList, "sessileburgspec", param->sessileburgspec, V_DBL, 30,
          VFLAG_NULL);

  BindVar(CPList, "sessilelinespec", param->sessilelinespec, V_DBL, 30,
          VFLAG_NULL);

  /*
   *      Flux decomposition
   */
  BindVar(CPList, "Flux decomposition", (void *)NULL, V_COMMENT, 1, VFLAG_NULL);

  BindVar(CPList, "totstraintensor", param->totstraintensor, V_DBL, 6,
          VFLAG_NULL);

  BindVar(CPList, "totpStn", param->totpStn, V_DBL, 6, VFLAG_NULL);

  BindVar(CPList, "totpSpn", param->totpSpn, V_DBL, 6, VFLAG_NULL);

  /*
   *      Total system density: this is informational only and recalculated
   *      each time a restart file is written.
   */
  BindVar(CPList, "Total density. Informational only; ignored on input",
          (void *)NULL, V_COMMENT, 1, VFLAG_NULL);

  BindVar(CPList, "disloDensity", &param->disloDensity, V_DBL, 1, VFLAG_NULL);

  /*
   *      Velocity statistics
   */
  BindVar(CPList, "Velocity statistics", (void *)NULL, V_COMMENT, 1,
          VFLAG_NULL);

  BindVar(CPList, "vAverage", &param->vAverage, V_DBL, 1, VFLAG_NULL);

  BindVar(CPList, "vStDev", &param->vStDev, V_DBL, 1, VFLAG_NULL);

  /*
   *      I/O controls and options
   */
  BindVar(CPList, "I/O controls and parameters", (void *)NULL, V_COMMENT, 1,
          VFLAG_NULL);

  BindVar(CPList, "dirname", param->dirname, V_STRING, 1, VFLAG_NULL);
  BindVar(CPList, "PrecipitatesFileName", param->PrecipitatesFileName, V_STRING,
          1, VFLAG_NULL);

  BindVar(CPList, "ShearablePrecipitates", &param->ShearablePrecipitates, V_INT,
          1, VFLAG_NULL);
  param->ShearablePrecipitates = 1;

  BindVar(CPList, "PairPrecipitateShearing", &param->PairPrecipitateShearing,
          V_INT, 1, VFLAG_NULL);
  param->PairPrecipitateShearing = 0;

  BindVar(CPList, "PrecipitateNucleationRate",
          &param->PrecipitateNucleationRate, V_DBL, 1, VFLAG_NULL);
  param->PrecipitateNucleationRate = 0.0;

  BindVar(CPList, "PrecipitateNucleationRadius",
          &param->PrecipitateNucleationRadius, V_DBL, 1, VFLAG_NULL);
  param->PrecipitateNucleationRadius = 0.0;

  BindVar(CPList, "PrecipitateNucleationMaxCount",
          &param->PrecipitateNucleationMaxCount, V_INT, 1, VFLAG_NULL);
  param->PrecipitateNucleationMaxCount = -1;

  BindVar(CPList, "APBFileName", param->APBFileName, V_STRING, 1, VFLAG_NULL);
  strcpy(param->APBFileName, "");

  BindVar(CPList, "numIOGroups", &param->numIOGroups, V_INT, 1, VFLAG_NULL);
  param->numIOGroups = 1;

  /*
   *      restart files
   */

  BindVar(CPList, "savefreq", &param->savefreq, V_INT, 1, VFLAG_NULL);
  param->savefreq = 100;

  BindVar(CPList, "savecounter", &param->savecounter, V_INT, 1, VFLAG_NULL);
  param->savecounter = 0;

  /*
   *      properties files
   */

  BindVar(CPList, "savepropdt", &param->savepropdt, V_DBL, 1, VFLAG_NULL);
  param->savepropdt = -1.0;

  BindVar(CPList, "saveproptime", &param->saveproptime, V_DBL, 1, VFLAG_NULL);

  /*
   *      tecplot files
   */
  BindVar(CPList, "tecplot", &param->tecplot, V_INT, 1, VFLAG_NULL);

  BindVar(CPList, "tecplotdt", &param->tecplotdt, V_DBL, 1, VFLAG_NULL);
  param->tecplotdt = -1.0;

  BindVar(CPList, "tecplottime", &param->tecplottime, V_DBL, 1, VFLAG_NULL);

  BindVar(CPList, "tecplotcounter", &param->tecplotcounter, V_INT, 1,
          VFLAG_NULL);

  /*
   *      3D-dislocation density data
   */
  BindVar(CPList, "savedensityspec", param->savedensityspec, V_INT, 3,
          VFLAG_NULL);

  /*
   *      Miscellaneous parameters
   */
  BindVar(CPList, "Miscellaneous parameters", (void *)NULL, V_COMMENT, 1,
          VFLAG_NULL);
}

void DataParamInit(Param_t *param, ParamList_t *DPList) {
  DPList->paramCnt = 0;
  DPList->varList = NULL;
  BindVar(DPList, "dataFileVersion", &param->dataFileVersion, V_INT, 1,
          VFLAG_NULL);
  param->dataFileVersion = NODEDATA_FILE_VERSION;
  // BindVar(DPList, "Center", param->center, V_DBL, 3,VFLAG_NULL);
  BindVar(DPList, "Dimensions", param->Dimensions, V_DBL, 3, VFLAG_NULL);
  /*cout <<"minX" <<param->Dimensions[0] <<endl;
  cout <<"maxX" <<param->Dimensions[3] <<endl;
  cout <<"minY" <<param->Dimensions[1] <<endl;
  cout <<"maxY" <<param->Dimensions[4] <<endl;
  cout <<"minZ" <<param->Dimensions[2] <<endl;
  cout <<"maxZ" <<param->Dimensions[5] <<endl;
  cout <<"Xlength" <<"\t" <<param->Dimensions[0] <<endl;
  cout <<"Ylength" <<"\t" <<param->Dimensions[1] <<endl;
  cout <<"Zlength" <<"\t" <<param->Dimensions[2] <<endl;
  param->Dimensions[0] = param->Dimensions[0]-param->center[0];
  param->Dimensions[1] = param->Dimensions[1]-param->center[1];
  param->Dimensions[2] = param->Dimensions[2]-param->center[2];
  param->Dimensions[3] = param->Dimensions[3]-param->center[0];
  param->Dimensions[4] = param->Dimensions[4]-param->center[1];
  param->Dimensions[5] = param->Dimensions[5]-param->center[2];
  param->Dimensions[0]= param->Dimensions[3]- param->Dimensions[0];
  param->Dimensions[1]= param->Dimensions[4]- param->Dimensions[1];
  param->Dimensions[2]= param->Dimensions[5]- param->Dimensions[2];*/
  BindVar(DPList, "XAxis", param->XAxis, V_DBL, 3, VFLAG_NULL);
  param->XAxis[0] = 1.0;
  param->XAxis[1] = 0.0;
  param->XAxis[2] = 0.0;
  BindVar(DPList, "YAxis", param->YAxis, V_DBL, 3, VFLAG_NULL);
  param->YAxis[0] = 0.0;
  param->YAxis[1] = 1.0;
  param->YAxis[2] = 0.0;
  BindVar(DPList, "nodeCount", &param->nodeCount, V_INT, 1, VFLAG_NULL);
  BindVar(DPList, "GrainDiameter", &param->GrainDiameter, V_DBL, 1, VFLAG_NULL);
  param->GrainDiameter = 4000.0;
  BindVar(DPList, "GeometryType", &param->GeometryType, V_INT, 1, VFLAG_NULL);
  param->GeometryType = 1;
}

void StressParamInit(Param_t *param, ParamList_t *SPList) {
  /*
   *      Note: Parameters need only be initialized if their
   *      default values are non-zero.
   */
  SPList->paramCnt = 0;
  SPList->varList = NULL;
  BindVar(SPList, "Time_evo", &param->Time_evo, V_INT, 1, VFLAG_NULL);
  param->Time_evo = 1;
  BindVar(SPList, "GP_x", &param->GP_x, V_INT, 1, VFLAG_NULL);
  param->GP_x = 2;
  BindVar(SPList, "GP_y", &param->GP_y, V_INT, 1, VFLAG_NULL);
  param->GP_y = 2;
  BindVar(SPList, "GP_z", &param->GP_z, V_INT, 1, VFLAG_NULL);
  param->GP_z = 2;
  BindVar(SPList, "stress_Dim", param->stress_Dim, V_DBL, 3, VFLAG_NULL);
}
