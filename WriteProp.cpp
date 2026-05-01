/****************************************************************************
 *
 *      Function:     WriteProp
 *      Description:  Depending on the input parameter, write a particular
 *                    time-dependent property out to that property's
 *                    diagnostic file, along with a timestamp.
 *
 *      NOTE: The contents of the file written below containing the
 *      dislocation density have changed over time.  The contents of
 *      the various versions are defined below.
 *
 *      Version 0 contents:
 *        - strain
 *        - dislocation density
 *        - deleted segment length
 *        - average velocity
 *        - std deviation of dislocation velocities
 *
 *      Version 1 contents: After the velocity std deviation, the
 *      following items were added:
 *        - version number
 *        - a set of values indicating the dislocation density for
 *          segments of specific groupings of burgers vectors.  These
 *          groupings differ for BCC and FCC mobility.
 *
 *          BCC groupings:
 *
 *          group #     burgers vector types
 *            0         [ 1 1 1] [-1-1-1]
 *            1         [-1 1 1] [ 1-1-1]
 *            2         [ 1-1 1] [-1 1-1]
 *            3         [ 1 1-1] [-1-1 1]
 *            4         [ 1 0 0] [-1 0 0]
 *                      [ 0 1 0] [ 0-1 0]
 *                      [ 0 0 1] [ 0 0-1]
 *
 *          FCC groupings:
 *
 *          group #     burgers vector types
 *            0         [ 1 1 0] [-1-1 0]
 *            1         [-1 1 0] [ 1-1 0]
 *            2         [ 1 0 1] [-1 0-1]
 *            3         [-1 0 1] [ 1 0-1]
 *            4         [ 0 1 1] [ 0-1-1]
 *            5         [ 0-1 1] [ 0 1-1]
 *            6         all others
 *
 *      Version 2 contents:
 *        - plastic strain added as first column of data
 *
 ****************************************************************************/

#include "Home.h"
#include "WriteProp.h"
#include "Util.h"
#include "Mobility.h"

#define DENSITY_FILE_VERSION 2

void WriteProp(Home_t *home, int property) {
  int i, numItems, tmpOffset, numBurgVectors;
  real8 al, am, an, amag, sigijk, epsijk, pstnijk, dpstnijk;
  real8 *localDensityVals = NULL;
  real8 *globalDensityVals = NULL;
  real8 totDensityChange[14];
  char fileName[256];
  FILE *fp;
  Param_t *param;

  param = home->param;

  al = param->edotdir[0];
  am = param->edotdir[1];
  an = param->edotdir[2];

  amag = sqrt(al * al + am * am + an * an);

  al /= amag;
  am /= amag;
  an /= amag;

  switch (property) {
  case DENSITY:

#ifdef PARALLEL
    numItems = param->numBurgGroups + 1;
    tmpOffset = 1;

    localDensityVals = (real8 *)malloc((numItems) * sizeof(real8));
    globalDensityVals = (real8 *)malloc((numItems) * sizeof(real8));

    for (i = 0; i < numItems; i++) {
      globalDensityVals[i] = 0.0;
    }

    localDensityVals[0] = param->delSegLength;

    for (i = 0; i < param->numBurgGroups; i++) {
      localDensityVals[i + tmpOffset] = 0.0;
    }

    MPI_Allreduce(localDensityVals, globalDensityVals, numItems, MPI_DOUBLE,
                  MPI_SUM, MPI_COMM_WORLD);
    param->delSegLength = globalDensityVals[0];

    free(localDensityVals);
    free(globalDensityVals);
#endif

    /*
     *          convert to density and write results to the file
     */
    param->delSegLength *= param->burgVolFactor;

    if (home->myDomain == 0) {
      snprintf(fileName, sizeof(fileName), "%s/density", DIR_PROPERTIES);
      fp = fopen(fileName, "a");

      pstnijk =
          param->totpStn[0] * al * al + param->totpStn[1] * am * am +
          param->totpStn[2] * an * an + 2.0 * param->totpStn[3] * am * an +
          2.0 * param->totpStn[4] * an * al + 2.0 * param->totpStn[5] * al * am;

      /*
       *              First print the standard stuff that's common
       *              regardless of which mobility (BCC, FCC, etc)
       *              is in use.
       */
      fprintf(fp, "%e %e %e %e %e %e\n", pstnijk, param->eRate * param->timeNow,
              param->disloDensity, param->delSegLength, param->vAverage,
              param->vStDev);
      fclose(fp);

      if (home->poPrecipitateServer != (ParadisPrecipitateServer *)NULL) {
        snprintf(fileName, sizeof(fileName), "%s/precipitates", DIR_PROPERTIES);
        bool bWriteHeader = false;
        FILE *fpCheck = fopen(fileName, "r");
        if (fpCheck == NULL) {
          bWriteHeader = true;
        } else {
          fseek(fpCheck, 0, SEEK_END);
          if (ftell(fpCheck) <= 0) {
            bWriteHeader = true;
          }
          fclose(fpCheck);
        }
        fp = fopen(fileName, "a");
        if (bWriteHeader) {
          fprintf(fp,
                  "# cycle timeNow precipitatesCount totalNucleatedCount "
                  "nucleationRate nucleationRadius nucleationMaxCount\n");
        }
        fprintf(fp, "%d %e %u %u %e %e %d\n", home->cycle, param->timeNow,
                home->poPrecipitateServer->GetPrecipitatesCount(),
                home->poPrecipitateServer->GetTotalNucleatedCount(),
                param->PrecipitateNucleationRate,
                param->PrecipitateNucleationRadius,
                param->PrecipitateNucleationMaxCount);
        fclose(fp);
      }
    }

    /*
     *          Reinitialize accumulated length of deleted segments
     *          after writing
     */
    param->delSegLength = 0.0;
    break;

  case EPS:
    if (home->myDomain == 0) {
      sigijk = param->appliedStress[0] * al * al +
               param->appliedStress[1] * am * am +
               param->appliedStress[2] * an * an +
               2.0 * param->appliedStress[3] * am * an +
               2.0 * param->appliedStress[4] * an * al +
               2.0 * param->appliedStress[5] * al * am;

      epsijk = param->totstraintensor[0] * al * al +
               param->totstraintensor[1] * am * am +
               param->totstraintensor[2] * an * an +
               2.0 * param->totstraintensor[3] * am * an +
               2.0 * param->totstraintensor[4] * an * al +
               2.0 * param->totstraintensor[5] * al * am;

      pstnijk =
          param->totpStn[0] * al * al + param->totpStn[1] * am * am +
          param->totpStn[2] * an * an + 2.0 * param->totpStn[3] * am * an +
          2.0 * param->totpStn[4] * an * al + 2.0 * param->totpStn[5] * al * am;

      snprintf(fileName, sizeof(fileName), "%s/time_Plastic_strain",
               DIR_PROPERTIES);
      fp = fopen(fileName, "a");
      fprintf(fp, "%e %e\n", param->timeNow, pstnijk);
      fclose(fp);

      if ((param->loadType == 1) || (param->loadType == 4) ||
          (param->loadType == 7) || (param->loadType == 8)) {
        snprintf(fileName, sizeof(fileName), "%s/stress_Plastic_strain",
                 DIR_PROPERTIES);
        fp = fopen(fileName, "a");
        fprintf(fp, "%e %e\n", pstnijk, sigijk);
        fclose(fp);

        snprintf(fileName, sizeof(fileName), "%s/stress_Total_strain",
                 DIR_PROPERTIES);
        fp = fopen(fileName, "a");
        switch (param->loadType) {
        case 1:
          fprintf(fp, "%e %e\n", epsijk, sigijk);
          break;
        case 4:
          fprintf(fp, "%e %e %e %d\n", param->netCyclicStrain, sigijk,
                  param->timeNow, param->numLoadCycle);
          break;
        case 7:
          fprintf(fp, "%e %e\n", epsijk, sigijk);
          break;
        case 8:
          fprintf(fp, "%e %e\n", epsijk, sigijk);
          break;
        }
        fclose(fp);
      }
    }
    break;

  case ALL_EPS:
    if (home->myDomain == 0) {
      snprintf(fileName, sizeof(fileName), "%s/alleps", DIR_PROPERTIES);
      fp = fopen(fileName, "a");

      fprintf(fp, "%d %e %e %e %e %e %e %e %e\n", home->cycle, param->timeNow,
              param->totpStn[0], param->totpStn[1], param->totpStn[2],
              param->totpStn[3], param->totpStn[4], param->totpStn[5],
              param->disloDensity);
      fclose(fp);
    }
    break;
  /*case CrossSlip:
  if (home->myDomain == 0) {
          snprintf(fileName, sizeof(fileName), "%s/CrossSlip",
                   DIR_PROPERTIES);
          fp=fopen(fileName, "a");

          fprintf(fp,"%d %e %e %e %e %e\n",
                  home->cycle, param->timeNow,
                  param->GlideShear,param->CrossSlipShear,param->EscaigGlide,
                  param->EsacigCrossSlip);
          fclose(fp);
      }
      break;*/

  case EPSDOT:
    if (home->myDomain == 0) {
      if (param->deltaTT > 1.e-20) {
        dpstnijk = param->delpStrain[0] * al * al +
                   param->delpStrain[1] * am * am +
                   param->delpStrain[2] * an * an +
                   2.0 * param->delpStrain[3] * am * an +
                   2.0 * param->delpStrain[4] * an * al +
                   2.0 * param->delpStrain[5] * al * am;

        snprintf(fileName, sizeof(fileName), "%s/epsdot", DIR_PROPERTIES);
        fp = fopen(fileName, "a");
        fprintf(fp, "%e %e\n", param->timeNow, fabs(dpstnijk / param->deltaTT));
        fclose(fp);
      }
    }
    break;

  case DENSITY_DELTA:
    /*
     *          Currently only supported for BCC mobilities
     */
    if (param->materialType == MAT_TYPE_BCC) {

      /*
       *              We track 2 values (density gain and density loss) for
       *              each of 7 burgers vectors
       */
      numBurgVectors = 7;
      numItems = numBurgVectors * 2;
#ifdef PARALLEL

      for (i = 0; i < numItems; i++) {
        totDensityChange[i] = 0.0;
      }

      MPI_Reduce(param->densityChange, totDensityChange, numItems, MPI_DOUBLE,
                 MPI_SUM, 0, MPI_COMM_WORLD);
#else
      for (i = 0; i < numItems; i++) {
        totDensityChange[i] = param->densityChange[i];
      }
#endif
      if (home->myDomain == 0) {
        real8 totGain, totLoss;

        for (i = 0; i < numBurgVectors; i++) {
          totGain += totDensityChange[i];
          totLoss += totDensityChange[numBurgVectors + i];
        }

        snprintf(fileName, sizeof(fileName), "%s/density_delta",
                 DIR_PROPERTIES);

        fp = fopen(fileName, "a");

        fprintf(fp, "%e %e %e ", param->eRate * param->timeNow, totGain,
                totLoss);

        for (i = 0; i < numItems; i++) {
          fprintf(fp, "%e ", totDensityChange[i]);
        }
        fprintf(fp, "\n");

        fclose(fp);
      }
      /*
       *              Zero out the accumulated values since we've written the
       *              data to disk.
       */
      for (i = 0; i < numItems; i++) {
        param->densityChange[i] = 0.0;
      }
    }
    break;

  default:
    Fatal("WriteProp: input parameter (property=%d) invalid", property);
    break;

  } /* switch() */

  return;
}
