/**************************************************************************
 *
 *      Author:  Moono Rhee
 *      Function: LoadCurve
 *
 *      Description: This subroutine defines the type of load curves.
 *                   Works only with the conventional x-y-z (global)
 *                   coordinate system.  If loading axis rotated, the
 *                   loading axis can be rotated, or One can rotate
 *                   the initial dislocation configuration to a
 *                   "laboratory" coordinate system.
 *
 *                   Types of load curves:
 *                      0  Creep
 *                      1  Constant strain test
 *                      2  Displacement-controlled
 *                      3  Junction unzipping jump test
 *                      4  Total strain controlled cyclic load
 *                      5  Plastic strain controlled cyclic load
 *                      6  Load-time curve
 *
 *      Last Modified:  01/03/2001 - original version
 *                      03/13/2003 - M. Rhee Removed anisotropic elastic
 *                                   constants.  Modified to include
 *                                   isotropic Hooke's law for arbitray
 *                                   loading.
 *                      11/11/2003 - MasatoH Implementation of loading axis
 *                                   rotation due to accumuration of
 *                                   material spin.  Instead of crystal
 *                                   system, lab frame is rotated in opposite
 *                                   way
 *                      06/23/2004 - M.Rhee Added strain decomposition and
 *                                   density flux decompostion.  Modified
 *                                   message passing calls for all decomposed
 *                                   strain/density info
 *                      07/12/2004 - Masato Strain contolled cyclic load
 *                                   is implemented.
 *
 ***************************************************************************/
#include "Home.h"
#include "Util.h"
#include <stdio.h>
#include <math.h>

/*
 *      Ss(): Sine for small angle i.e. Ss~x
 *      Cs(): Cosine for small angle i.e. Cs~1-x^2/2
 */
#define Ss(a) ((a))
#define Cs(a) (1.0 - (0.5 * (a) * (a)))

/*
 *      Function:     SpinMatrix
 *      Description:  Small rotation matrix for accumulaed
 *                    rotations around axis 1, 2 and 3.
 *
 *                    Cs() = Cosine for small angle
 *                    Ss() = Sin for small angle
 */
static void SpinMatrix(real8 p1, real8 p2, real8 p3, real8 Rspin[3][3]) {
  Rspin[0][0] = Cs(p3) * Cs(p2);
  Rspin[1][1] = Cs(p3) * Cs(p1) + Ss(p3) * Ss(p2) * Ss(p1);
  Rspin[2][2] = Cs(p2) * Cs(p1);
  Rspin[0][1] = -Ss(p3) * Cs(p1) + Cs(p3) * Ss(p1) * Ss(p2);
  Rspin[1][2] = -Cs(p3) * Ss(p1) + Ss(p3) * Ss(p2) * Cs(p1);
  Rspin[2][0] = -Ss(p2);
  Rspin[0][2] = Ss(p3) * Ss(p1) + Cs(p3) * Cs(p1) * Ss(p2);
  Rspin[1][0] = Ss(p3) * Cs(p2);
  Rspin[2][1] = Cs(p2) * Ss(p1);
}

void LoadCurve(Home_t *home) {
  int i, j, k, loadtype, indxerate;
  int numLoadCycle, numLoadCycle2;
  real8 youngs, erate, dtt;
  real8 shr;
  real8 modulus, dpl_stn, dStress, amag, al, am, an;
  real8 sigijk, stn_cut;
  real8 phi1, phi2, phi3;
  real8 Rspin[3][3];
  real8 tempedot[3], temppassedot[3];
  real8 pstnijk, eAmp, timeNow, cTimeOld;
  real8 dCyclicStrain;
  real8 totCyclicStrain, netCyclicStrain;
  Param_t *param;

  param = home->param;
  loadtype = param->loadType;
  shr = param->shearModulus;
  youngs = 2.0 * shr * (1.0 + param->pois);
  erate = param->eRate;
  dtt = param->deltaTT;
  indxerate = param->indxErate;

  sigijk = 0;
  totCyclicStrain = 0.0;

  /*
   *      for cyclic load
   */
  eAmp = param->eAmp;
  timeNow = param->timeNow;
  cTimeOld = param->cTimeOld;
  numLoadCycle = param->numLoadCycle;
  netCyclicStrain = param->netCyclicStrain;
  dCyclicStrain = param->dCyclicStrain;

  /*
   *      If we're including osmotic forces on dislocation segments, we
   *      need to use the delta plastic strain to adjust the vacancy
   *      concentration.  Basically, sum the diagonal of the delta plastic
   *      strain and add to the vacancy concentration.
   */
  if (param->vacancyConcEquilibrium > 0.0) {
    param->vacancyConc +=
        (param->delpStrain[0] + param->delpStrain[1] + param->delpStrain[2]);
  }

  /*
   *      Part for Loading Axis Rotation due to Small Deformation Spin.
   *
   *      Some changes in rotation angles due to the deformation spins
   *      around x, y, and z axis
   */
  phi1 = -param->delpSpin[3];
  phi2 = param->delpSpin[4];
  phi3 = -param->delpSpin[5];

  /*
   *      Matrix for (combined) rotation around x,y,and z in the sequence.
   *      This sequential rotation is correct only for small changes
   *      in the angles since real material rotation occurs simultaneously.
   *      For counter-rotation, sign of phi is flipped.
   */
  SpinMatrix(phi1, phi2, phi3, Rspin);

  /*
   *      Compute nodal velocity : bug is fixed. Vector is address
   */
  tempedot[0] = param->edotdir[0];
  tempedot[1] = param->edotdir[1];
  tempedot[2] = param->edotdir[2];

  temppassedot[0] = 0.0;
  temppassedot[1] = 0.0;
  temppassedot[2] = 0.0;

  Matrix33Vector3Multiply(Rspin, tempedot, temppassedot);

  param->edotdir[0] = temppassedot[0];
  param->edotdir[1] = temppassedot[1];
  param->edotdir[2] = temppassedot[2];

  /*
   *      Arbitrary loading direction but keep in lab frame
   */
  al = param->edotdir[0];
  am = param->edotdir[1];
  an = param->edotdir[2];

  amag = sqrt(al * al + am * am + an * an);

  al /= amag;
  am /= amag;
  an /= amag;

  double dAmplitude = param->CyclicLoadAmplitude;
  double dTotalStrainIncrement = 0.0;
  switch (loadtype) {
  /*
   *          creep - what we have been using so this should be default
   */
  case 0:
    break;
  /*
   *          constant strain rate
   */
  case 1: {
    /*
     *              Cover for specific loading direction also
     */
    dpl_stn =
        param->delpStrain[0] * al * al + param->delpStrain[1] * am * am +
        param->delpStrain[2] * an * an +
        2.0 * (param->delpStrain[3] * am * an + param->delpStrain[4] * an * al +
               param->delpStrain[5] * al * am);

    if (indxerate <= 3) {
      modulus = youngs;
    } else {
      modulus = 2.0 * shr;
    }

    /*
     *              local in the [l m n] frame
     */
    dStress = modulus * (erate * dtt - dpl_stn);
    /*
     *              global (100)-(010)-(001) frame
     */
    param->appliedStress[0] += dStress * al * al;
    param->appliedStress[1] += dStress * am * am;
    param->appliedStress[2] += dStress * an * an;
    param->appliedStress[3] += dStress * an * am;
    param->appliedStress[4] += dStress * an * al;
    param->appliedStress[5] += dStress * al * am;

    param->totstraintensor[0] += erate * dtt * al * al;
    param->totstraintensor[1] += erate * dtt * am * am;
    param->totstraintensor[2] += erate * dtt * an * an;
    param->totstraintensor[3] += erate * dtt * an * am;
    param->totstraintensor[4] += erate * dtt * an * al;
    param->totstraintensor[5] += erate * dtt * al * am;

    break;
  }
  /*
   *          jump test
   */
  case 2: {
    //                 stn_cut = 5.e-12;
    //                 dStress = 1.e5;
    //                 dpl_stn = param->delpStrain[indxerate-1];
    //
    //                 if ((dpl_stn > 0.0) && (dpl_stn < stn_cut))
    //                 {
    //                     param->appliedStress[indxerate-1] += dStress;
    //                     home->cycle = 1;
    //                 }
    //                 else if (dpl_stn < 0)
    //                 {
    //
    //                 }
    //                 else
    //                 {
    //                     if (home->cycle > 200000)
    //                     {
    //                          printf("Must be ok now printing the critical
    //                          stress \n"); printf("Critical Stress =%e \n",
    //                          sigijk); printf("  minSeg   =%e \n",
    //                          param->minSeg); printf("  maxSeg   =%e \n",
    //                          param->maxSeg); Fatal("Doing clean terminate in
    //                          LoadCurve");
    //                     }
    //                }
    //
    //                home->cycle++;
    //
    //                sigijk  =  param->appliedStress[indxerate-1];
    //                dpl_stn =  param->delpStrain[indxerate-1];
    //
    //                if (home->cycle % 50 == 0)
    //                {
    //                    printf("sig=%e stn=%e dt=%e cyc=%d\n", sigijk,dpl_stn,
    //                    param->deltaTT, home->cycle);
    //                }

    break;
  }
  /*
   *          Junction unzipping jump test; not for general case yet */
  case 3: {
#if 1
    stn_cut = 1.0e-10; /* lowbound strain cut for numerical */
                       /* stability                         */
    dStress = 1.e5;
    dpl_stn = param->delpStrain[indxerate - 1];

    if ((dpl_stn > 0.0) && (dpl_stn < stn_cut)) {
      param->appliedStress[indxerate - 1] += dStress;
    }

    sigijk = param->appliedStress[indxerate - 1];
    dpl_stn = param->delpStrain[indxerate - 1];

    if (home->cycle % 10 == 0) {
      printf("sig=%e stn=%e dt=%e cyc=%d\n", sigijk, dpl_stn, param->deltaTT,
             home->cycle);
    }
#endif
    break;
  }

  /*
   *          strain control cyclic load
   *
   *              stainCycle    = current loading cycle
   *              eAmp          = strain amplitude for each side
   *              dCyclicStrain = change in the strain for each side
   *              acumStrain    = accumulated strain
   *              sgnLoad       = sign of load
   */
  case 4: {

    /*
     *              Cover for specific loading direction also
     */
    dpl_stn = param->delpStrain[0] * al * al + param->delpStrain[1] * am * am +
              param->delpStrain[2] * an * an +
              2.0 * param->delpStrain[3] * am * an +
              2.0 * param->delpStrain[4] * an * al +
              2.0 * param->delpStrain[5] * al * am;

    if (indxerate <= 3) {
      modulus = youngs;
    } else {
      modulus = 2.0 * shr;
    }

    dCyclicStrain = erate * dtt;
    param->dCyclicStrain = dCyclicStrain;

    /*
     *              local in the [l m n] frame
     */
    dStress = modulus * (dCyclicStrain - dpl_stn);

    totCyclicStrain = fabs(erate * timeNow);
    numLoadCycle = rint(0.5 * totCyclicStrain / eAmp);
    numLoadCycle2 = rint(0.5 * totCyclicStrain / eAmp - 0.5);

    netCyclicStrain = fmod(totCyclicStrain, 2 * eAmp);

    if (fabs(netCyclicStrain) > eAmp) {
      netCyclicStrain = 2 * eAmp - netCyclicStrain;
    }

    netCyclicStrain = pow(-1.0, numLoadCycle2) * fabs(netCyclicStrain);

    param->netCyclicStrain = netCyclicStrain;

    cTimeOld = timeNow;
    erate = fabs(erate) * pow(-1.0, numLoadCycle);
    dCyclicStrain = 0;
    param->cTimeOld = cTimeOld;
    param->numLoadCycle = numLoadCycle;
    param->eRate = erate;

    /*
     *              global (100)-(010)-(001) frame
     */
    param->appliedStress[0] += dStress * al * al;
    param->appliedStress[1] += dStress * am * am;
    param->appliedStress[2] += dStress * an * an;
    param->appliedStress[3] += dStress * an * am;
    param->appliedStress[4] += dStress * an * al;
    param->appliedStress[5] += dStress * al * am;

    /*
                    param->totstraintensor[0] = erate * param->timeNow *  al*al;
                    param->totstraintensor[1] = erate * param->timeNow *  am*am;
                    param->totstraintensor[2] = erate * param->timeNow *  an*an;
                    param->totstraintensor[3] = erate * param->timeNow *  an*am;
                    param->totstraintensor[4] = erate * param->timeNow *  an*al;
                    param->totstraintensor[5] = erate * param->timeNow *  al*am;
    */
    param->totstraintensor[0] = netCyclicStrain * al * al;

    break;
  }

  /*
   *             Plastic strain control cyclic load
   *                 stainCycle    = current loading cycle
   *                 eAmp          = strain amplitude for each side
   *                 dCyclicStrain = change in the strain for each side
   *                 acumStrain    = accumulated strain
   *                 sgnLoad       = sign of load
   */
  case 5: {

    /*
     *              Cover for specific loading direction als0
     */
    pstnijk = param->totpStn[0] * al * al + param->totpStn[1] * am * am +
              param->totpStn[2] * an * an + 2.0 * param->totpStn[3] * am * an +
              2.0 * param->totpStn[4] * an * al +
              2.0 * param->totpStn[5] * al * am;

    dpl_stn = param->delpStrain[0] * al * al + param->delpStrain[1] * am * am +
              param->delpStrain[2] * an * an +
              2.0 * param->delpStrain[3] * am * an +
              2.0 * param->delpStrain[4] * an * al +
              2.0 * param->delpStrain[5] * al * am;

    if (indxerate <= 3) {
      modulus = youngs;
    } else {
      modulus = 2.0 * shr;
    }

    dCyclicStrain = erate * dtt;
    param->dCyclicStrain = dCyclicStrain;

    /*
     *              local in the [l m n] frame
     */
    dStress = modulus * (dCyclicStrain - dpl_stn);

    totCyclicStrain += fabs(dpl_stn);
    numLoadCycle = rint(0.5 * pstnijk / eAmp);
    numLoadCycle2 = rint(0.5 * pstnijk / eAmp - 0.5);

    netCyclicStrain = fmod(pstnijk, 2 * eAmp);

    if (fabs(netCyclicStrain) > eAmp) {
      netCyclicStrain = 2 * eAmp - netCyclicStrain;
    }

    netCyclicStrain = pow(-1.0, numLoadCycle2) * fabs(netCyclicStrain);

    param->netCyclicStrain = netCyclicStrain;

    cTimeOld = timeNow;
    erate = fabs(erate) * pow(-1.0, numLoadCycle);
    dCyclicStrain = 0;
    param->cTimeOld = cTimeOld;
    param->numLoadCycle = numLoadCycle;
    param->eRate = erate;

    /*
     *              global (100)-(010)-(001) frame
     */
    param->appliedStress[0] += dStress * al * al;
    param->appliedStress[1] += dStress * am * am;
    param->appliedStress[2] += dStress * an * an;
    param->appliedStress[3] += dStress * an * am;
    param->appliedStress[4] += dStress * an * al;
    param->appliedStress[5] += dStress * al * am;

    /*
                    param->totstraintensor[0] = erate * param->timeNow *  al*al;
                    param->totstraintensor[1] = erate * param->timeNow *  am*am;
                    param->totstraintensor[2] = erate * param->timeNow *  an*an;
                    param->totstraintensor[3] = erate * param->timeNow *  an*am;
                    param->totstraintensor[4] = erate * param->timeNow *  an*al;
                    param->totstraintensor[5] = erate * param->timeNow *  al*am;
    */
    param->totstraintensor[0] = netCyclicStrain * al * al;
    param->totstraintensor[1] = netCyclicStrain * am * am;
    param->totstraintensor[2] = netCyclicStrain * an * an;
    param->totstraintensor[3] = netCyclicStrain * an * am;
    param->totstraintensor[4] = netCyclicStrain * an * al;
    param->totstraintensor[5] = netCyclicStrain * al * am;

    break;
  }
  /*
   *          User defined load-time curve
   */
  case 6:
    break;

  case 7: {
    dpl_stn =
        param->delpStrain[0] * al * al + param->delpStrain[1] * am * am +
        param->delpStrain[2] * an * an +
        2.0 * (param->delpStrain[3] * am * an + param->delpStrain[4] * an * al +
               param->delpStrain[5] * al * am);
    modulus = youngs;
    double dT = 1.0 / param->CyclicLoadFrequency;
    double dT1 = param->timeNow / dT;
    double dT2 = dT1 - floor(dT1);
    double dSlope =
        4.0 * param->CyclicLoadAmplitude * param->CyclicLoadFrequency;
    if ((dT2 >= 0.25) && (dT2 < 0.75)) {
      dSlope = -dSlope;
    }

    // local in the [l m n] frame
    dTotalStrainIncrement = dSlope * dtt;
    dStress = modulus * (dTotalStrainIncrement - dpl_stn);

    // global (100)-(010)-(001) frame
    param->appliedStress[0] += dStress * al * al;
    param->appliedStress[1] += dStress * am * am;
    param->appliedStress[2] += dStress * an * an;
    param->appliedStress[3] += dStress * an * am;
    param->appliedStress[4] += dStress * an * al;
    param->appliedStress[5] += dStress * al * am;

    param->totstraintensor[0] += dTotalStrainIncrement * al * al;
    param->totstraintensor[1] += dTotalStrainIncrement * am * am;
    param->totstraintensor[2] += dTotalStrainIncrement * an * an;
    param->totstraintensor[3] += dTotalStrainIncrement * an * am;
    param->totstraintensor[4] += dTotalStrainIncrement * an * al;
    param->totstraintensor[5] += dTotalStrainIncrement * al * am;

    break;
  }

  case 8: {
    dpl_stn =
        param->delpStrain[0] * al * al + param->delpStrain[1] * am * am +
        param->delpStrain[2] * an * an +
        2.0 * (param->delpStrain[3] * am * an + param->delpStrain[4] * an * al +
               param->delpStrain[5] * al * am);
    modulus = youngs;

    /* local in the [l m n] frame */
    dTotalStrainIncrement = erate * dtt;
    dStress = modulus * (dTotalStrainIncrement - dpl_stn);
    // if the plastic strain in the current step is higher than the required
    // increment, don't increase the stress
    if (dTotalStrainIncrement < dpl_stn) {
      dStress = 0.0;
      dTotalStrainIncrement = dpl_stn;
    }

    double dIdealPreviousTotalStrain = erate * param->timeNow - erate * dtt;
    double dPreviousTotalStrain = param->totstraintensor[0] * al * al +
                                  param->totstraintensor[1] * am * am +
                                  param->totstraintensor[2] * an * an +
                                  2.0 * (param->totstraintensor[3] * am * an +
                                         param->totstraintensor[4] * an * al +
                                         param->totstraintensor[5] * al * am);
    if (dIdealPreviousTotalStrain < dPreviousTotalStrain) {
      dStress = 0.0;
      dTotalStrainIncrement = dpl_stn;
    }

    /* global (100)-(010)-(001) frame */
    param->appliedStress[0] += dStress * al * al;
    param->appliedStress[1] += dStress * am * am;
    param->appliedStress[2] += dStress * an * an;
    param->appliedStress[3] += dStress * an * am;
    param->appliedStress[4] += dStress * an * al;
    param->appliedStress[5] += dStress * al * am;

    param->totstraintensor[0] += dTotalStrainIncrement * al * al;
    param->totstraintensor[1] += dTotalStrainIncrement * am * am;
    param->totstraintensor[2] += dTotalStrainIncrement * an * an;
    param->totstraintensor[3] += dTotalStrainIncrement * an * am;
    param->totstraintensor[4] += dTotalStrainIncrement * an * al;
    param->totstraintensor[5] += dTotalStrainIncrement * al * am;

    break;
  }

  default: {
    // Fatal("Load curves not defined. Stopping the program. \n");
    break;
  }
  }
}
