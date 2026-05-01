/***************************************************************************
 * Function	: deWitInteraction
 * Description	: Calculates interaction stresses
 * 10/23/01 M.Rhee
 * 03/06/02 T.Pierce - modified the calling sequence, to pass in the
 *                     the coordinates, already periodically shifted as
 *                     necessary
 * 05/15/03 M.Hiratani - moved conditional statement out of Interpolation
 *                       removed unnecessary cal. for box size
 * 09/16/03 T.Pierce - Unrolled a bunch of loops in dSegImgStress for
 *                     optimization
 * 05/29/04 G.Hommes -	Removed obsolete functions calplanestress() and
 *			CalStress_fullimage().
 ***************************************************************************/
#include "Home.h"
#include "Util.h"
#include <math.h>

#ifdef PARALLEL
#include <mpi.h>
#endif

/* indices of 3rd derivatives */

/* PBC stress table */
#ifndef MAXSTRGRID
#define MAXSTRGRID 51
#endif

static int Offset1, Offset2, Offset3;
static real8 *RIJMTABLE;
static real8 *RIJMPBCTABLE;
static real8 RIJMLX, RIJMLY, RIJMLZ;

void FreeRijm() {
  if (RIJMTABLE != NULL) {
    free(RIJMTABLE);
    RIJMTABLE = NULL;
  }
}

void FreeRijmPBC() {
  if (RIJMPBCTABLE != NULL) {
    free(RIJMPBCTABLE);
    RIJMPBCTABLE = NULL;
  }
}

void ReadRijm(Home_t *home) {
  FILE *fp;
  Param_t *param;
  int i, j, k, m;
  int NX, NY, NZ;
  int NIMGX, NIMGY, NIMGZ, RIJM_Table_Size;
  real8 Lx, Ly, Lz;
  real8 ratiox, ratioy, ratioz;
  real8 rijml[3];

  param = home->param;

  Offset1 = 10 * MAXSTRGRID * MAXSTRGRID;
  Offset2 = 10 * MAXSTRGRID;
  Offset3 = 10;

  RIJM_Table_Size = (MAXSTRGRID * MAXSTRGRID * MAXSTRGRID * 10) * sizeof(real8);
  RIJMTABLE = (real8 *)calloc(1, RIJM_Table_Size);
  /*
   *	Only processor zero reads in the table.  It will broadcast the
   *	table (and supporting data) to all other processors
   */
  if (home->myDomain == 0) {
    printf("Reading file %s ...\n", param->Rijmfile);
    fp = fopen(param->Rijmfile, "r");
    if (fp == NULL) {
      fprintf(stderr, "ReadRijm file (%s) not found!\n", param->Rijmfile);
      exit(-1);
    }

    fscanf(fp, "%d %d %d", &NX, &NY, &NZ);
    fscanf(fp, "%le %le %le", &Lx, &Ly, &Lz);
    fscanf(fp, "%d %d %d", &NIMGX, &NIMGY, &NIMGZ);

    if ((NX > MAXSTRGRID) || (NY > MAXSTRGRID) || (NZ > MAXSTRGRID)) {
      fprintf(stderr, "NX=%d NY=%d NZ=%d MAXSTRGRID=%d\n", NX, NY, NZ,
              MAXSTRGRID);
      Fatal("table size larger than static array");
    }

    RIJMLX = Lx;
    RIJMLY = Ly;
    RIJMLZ = Lz;

    Lx = param->Lx;
    Ly = param->Ly;
    Lz = param->Lz;

    ratiox = RIJMLX / Lx;
    ratioy = RIJMLY / Ly;
    ratioz = RIJMLZ / Lz;

    if ((fabs(ratiox - ratioy) > 1e-3) || (fabs(ratioy - ratioz) > 1e-3)) {
      printf("Lx=%e Ly=%e Lz=%e\nLX=%e LY=%e LZ=%e\n", Lx, Ly, Lz, RIJMLX,
             RIJMLY, RIJMLZ);
      Fatal("Wrong RIJM Table read in");
    }

    for (i = 0; i < NX; i++) {
      for (j = 0; j < NY; j++) {
        for (k = 0; k < NZ; k++) {
          for (m = 0; m < 10; m++) {
            fscanf(fp, "%le\n",
                   &RIJMTABLE[i * Offset1 + j * Offset2 + k * Offset3 + m]);
          }
        }
      }
    }

    fclose(fp);

    param->imgstrgrid[0] = NX;
    param->imgstrgrid[1] = NY;
    param->imgstrgrid[2] = NZ;

    param->imgstrgrid[3] = NIMGX;
    param->imgstrgrid[4] = NIMGY;
    param->imgstrgrid[5] = NIMGZ;

    rijml[0] = RIJMLX;
    rijml[1] = RIJMLY;
    rijml[2] = RIJMLZ;

    printf("done.\n");
  }

#ifdef PARALLEL
  MPI_Bcast((char *)RIJMTABLE, RIJM_Table_Size, MPI_CHAR, 0, MPI_COMM_WORLD);
  MPI_Bcast((char *)&param->imgstrgrid[0], sizeof(param->imgstrgrid), MPI_CHAR,
            0, MPI_COMM_WORLD);
  MPI_Bcast(&rijml[0], 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif

  RIJMLX = rijml[0];
  RIJMLY = rijml[1];
  RIJMLZ = rijml[2];
}

void ReadRijmPBC(Home_t *home) {
  FILE *fp;
  Param_t *param;
  int i, j, k, m;
  int NX, NY, NZ;
  int NIMGX, NIMGY, NIMGZ, RIJM_Table_Size;
  real8 Lx, Ly, Lz;
  real8 ratiox, ratioy, ratioz;
  real8 rijml[3];

  param = home->param;

  Offset1 = 10 * MAXSTRGRID * MAXSTRGRID;
  Offset2 = 10 * MAXSTRGRID;
  Offset3 = 10;

  RIJM_Table_Size = (MAXSTRGRID * MAXSTRGRID * MAXSTRGRID * 10) * sizeof(real8);
  RIJMPBCTABLE = (real8 *)calloc(1, RIJM_Table_Size);
  /*
   *	Only processor zero reads in the table.  It will broadcast the
   *	table (and supporting data) to all other processors
   */
  if (home->myDomain == 0) {

    printf("Reading file %s ...\n", param->RijmPBCfile);
    fp = fopen(param->RijmPBCfile, "r");
    if (fp == NULL) {
      fprintf(stderr, "ReadRijmPBC file (%s) not found!\n", param->RijmPBCfile);
      Fatal("ReadRijmPBC() error");
    }

    fscanf(fp, "%d %d %d", &NX, &NY, &NZ);
    fscanf(fp, "%le %le %le", &Lx, &Ly, &Lz);
    fscanf(fp, "%d %d %d", &NIMGX, &NIMGY, &NIMGZ);

    RIJMLX = Lx;
    RIJMLY = Ly;
    RIJMLZ = Lz;

    Lx = param->Lx;
    Ly = param->Ly;
    Lz = param->Lz;

    ratiox = RIJMLX / Lx;
    ratioy = RIJMLY / Ly;
    ratioz = RIJMLZ / Lz;

    if ((fabs(ratiox - ratioy) > 1e-3) || (fabs(ratioy - ratioz) > 1e-3)) {
      printf("Lx=%e Ly=%e Lz=%e\n"
             "LX=%e LY=%e LZ=%e\n",
             Lx, Ly, Lz, RIJMLX, RIJMLY, RIJMLZ);
      Fatal("Wrong RIJMPBC Table read in");
    }

    for (i = 0; i < NX; i++)
      for (j = 0; j < NY; j++)
        for (k = 0; k < NZ; k++)
          for (m = 0; m < 10; m++) {
            fscanf(fp, "%le\n",
                   &RIJMPBCTABLE[i * Offset1 + j * Offset2 + k * Offset3 + m]);
          }

    fclose(fp);

    param->imgstrgrid[0] = NX;
    param->imgstrgrid[1] = NY;
    param->imgstrgrid[2] = NZ;

    param->imgstrgrid[3] = NIMGX;
    param->imgstrgrid[4] = NIMGY;
    param->imgstrgrid[5] = NIMGZ;

    rijml[0] = RIJMLX;
    rijml[1] = RIJMLY;
    rijml[2] = RIJMLZ;

    printf("done.\n");
  }

#ifdef PARALLEL
  MPI_Bcast((char *)RIJMPBCTABLE, RIJM_Table_Size, MPI_CHAR, 0, MPI_COMM_WORLD);
  MPI_Bcast((char *)&param->imgstrgrid[0], sizeof(param->imgstrgrid), MPI_CHAR,
            0, MPI_COMM_WORLD);
  MPI_Bcast(&rijml[0], 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif

  RIJMLX = rijml[0];
  RIJMLY = rijml[1];
  RIJMLZ = rijml[2];
}

void InterpolateRijm(Home_t *home, real8 rt[3], real8 Rijmarray[10], int pbc) {
  Param_t *param;
  real8 Lx, Ly, Lz, alpha, beta, gamma, xmin, ymin, zmin;
  real8 ratiox, ratio2;
  int m, NX, NY, NZ, nx, ny, nz;

  param = home->param;
  NX = param->imgstrgrid[0];
  NY = param->imgstrgrid[1];
  NZ = param->imgstrgrid[2];

  xmin = -0.5 * param->Dimensions[X];
  ymin = -0.5 * param->Dimensions[Y];
  zmin = -0.5 * param->Dimensions[Z];

  Lx = param->Lx;
  Ly = param->Ly;
  Lz = param->Lz;
  ratiox = RIJMLX / Lx;
  ratio2 = ratiox * ratiox;

  alpha = (rt[0] - xmin) / Lx * (NX - 1);
  beta = (rt[1] - ymin) / Ly * (NY - 1);
  gamma = (rt[2] - zmin) / Lz * (NZ - 1);

  while (alpha < 0)
    alpha += (NX - 1);
  while (alpha >= NX - 1)
    alpha -= (NX - 1);
  while (beta < 0)
    beta += (NY - 1);
  while (beta >= NY - 1)
    beta -= (NY - 1);
  while (gamma < 0)
    gamma += (NZ - 1);
  while (gamma >= NZ - 1)
    gamma -= (NZ - 1);

#ifdef _CUBICINTERPOLATE
  nx = (int)floor(alpha);
  ny = (int)floor(beta);
  nz = (int)floor(gamma);

  alpha -= nx;
  beta -= ny;
  gamma -= nz;

  if (pbc)
    for (m = 0; m < 10; m++)
      Rijmarray[m] =
          (1 - alpha) * (1 - beta) * (1 - gamma) *
              RIJMPBCTABLE[nx * Offset1 + ny * Offset2 + nz * Offset3 + m] +
          (1 - alpha) * (1 - beta) *
              (gamma)*RIJMPBCTABLE[nx * Offset1 + ny * Offset2 +
                                   (nz + 1) * Offset3 + m] +
          (1 - alpha) * (beta) * (1 - gamma) *
              RIJMPBCTABLE[nx * Offset1 + (ny + 1) * Offset2 + nz * Offset3 +
                           m] +
          (1 - alpha) * (beta) *
              (gamma)*RIJMPBCTABLE[nx * Offset1 + (ny + 1) * Offset2 +
                                   (nz + 1) * Offset3 + m] +
          (alpha) * (1 - beta) * (1 - gamma) *
              RIJMPBCTABLE[(nx + 1) * Offset1 + ny * Offset2 + nz * Offset3 +
                           m] +
          (alpha) * (1 - beta) *
              (gamma)*RIJMPBCTABLE[(nx + 1) * Offset1 + ny * Offset2 +
                                   (nz + 1) * Offset3 + m] +
          (alpha) * (beta) * (1 - gamma) *
              RIJMPBCTABLE[(nx + 1) * Offset1 + (ny + 1) * Offset2 +
                           nz * Offset3 + m] +
          (alpha) * (beta) *
              (gamma)*RIJMPBCTABLE[(nx + 1) * Offset1 + (ny + 1) * Offset2 +
                                   (nz + 1) * Offset3 + m];
  else
    for (m = 0; m < 10; m++)
      Rijmarray[m] =
          (1 - alpha) * (1 - beta) * (1 - gamma) *
              RIJMTABLE[nx * Offset1 + ny * Offset2 + nz * Offset3 + m] +
          (1 - alpha) * (1 - beta) *
              (gamma)*RIJMTABLE[nx * Offset1 + ny * Offset2 +
                                (nz + 1) * Offset3 + m] +
          (1 - alpha) * (beta) * (1 - gamma) *
              RIJMTABLE[nx * Offset1 + (ny + 1) * Offset2 + nz * Offset3 + m] +
          (1 - alpha) * (beta) *
              (gamma)*RIJMTABLE[nx * Offset1 + (ny + 1) * Offset2 +
                                (nz + 1) * Offset3 + m] +
          (alpha) * (1 - beta) * (1 - gamma) *
              RIJMTABLE[(nx + 1) * Offset1 + ny * Offset2 + nz * Offset3 + m] +
          (alpha) * (1 - beta) *
              (gamma)*RIJMTABLE[(nx + 1) * Offset1 + ny * Offset2 +
                                (nz + 1) * Offset3 + m] +
          (alpha) * (beta) * (1 - gamma) *
              RIJMTABLE[(nx + 1) * Offset1 + (ny + 1) * Offset2 + nz * Offset3 +
                        m] +
          (alpha) * (beta) *
              (gamma)*RIJMTABLE[(nx + 1) * Offset1 + (ny + 1) * Offset2 +
                                (nz + 1) * Offset3 + m];

#else
  nx = (int)rint(alpha);
  ny = (int)rint(beta);
  nz = (int)rint(gamma);

  if (pbc)
    for (m = 0; m < 10; m++)
      Rijmarray[m] =
          RIJMPBCTABLE[nx * Offset1 + ny * Offset2 + nz * Offset3 + m];
  else
    for (m = 0; m < 10; m++)
      Rijmarray[m] = RIJMTABLE[nx * Offset1 + ny * Offset2 + nz * Offset3 + m];

#endif
  /* scale stress */
  for (m = 0; m < 10; m++)
    Rijmarray[m] *= ratio2;
}

void dSegImgStress(Home_t *home, real8 Sigma[3][3], real8 px, real8 py,
                   real8 pz, real8 dlx, real8 dly, real8 dlz, real8 burgX,
                   real8 burgY, real8 burgZ, real8 rx, real8 ry, real8 rz,
                   int pbc) {
  ////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////
  Sigma[0][0] = 0.0;
  Sigma[0][1] = 0.0;
  Sigma[0][2] = 0.0;
  Sigma[1][0] = 0.0;
  Sigma[1][1] = 0.0;
  Sigma[1][2] = 0.0;
  Sigma[2][0] = 0.0;
  Sigma[2][1] = 0.0;
  Sigma[2][2] = 0.0;
  return;
  ////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////
  /* calculate image stress due to a differential dislocation segment
   * using the table Rijm
   * Sigma[3][3] return stress
   * (px, py, pz) location of segment (center)
   * (dlx, dly, dlz) vector of segment length
   * (burgX, burgY, burgZ) Burgers vector
   * (rx, ry, rz) location of field point
   *
   * Wei Cai, 7/6/2002
   */
  real8 rt[3], Rmpp[3], dl[3], bxdl[3], bxRmpp[3];
  real8 Rijmarray[10], Rijm[3][3][3], Rijmterm[3][3][3];
  real8 MUover4pi, onepois;
  int i, j, m;

  MUover4pi = home->param->shearModulus / 4 / M_PI;
  onepois = 1. / (1. - home->param->pois);

  dl[0] = dlx;
  dl[1] = dly;
  dl[2] = dlz;
  rt[0] = rx - px;
  rt[1] = ry - py;
  rt[2] = rz - pz;

  InterpolateRijm(home, rt, Rijmarray, pbc);

  /* unroll above loop */

  Rijm[0][0][0] = Rijmarray[0];
  Rijm[1][1][1] = Rijmarray[1];
  Rijm[2][2][2] = Rijmarray[2];
  Rijm[0][0][1] = Rijmarray[3];
  Rijm[0][1][0] = Rijmarray[3];
  Rijm[1][0][0] = Rijmarray[3];
  Rijm[0][0][2] = Rijmarray[4];
  Rijm[0][2][0] = Rijmarray[4];
  Rijm[2][0][0] = Rijmarray[4];
  Rijm[0][1][1] = Rijmarray[5];
  Rijm[1][0][1] = Rijmarray[5];
  Rijm[1][1][0] = Rijmarray[5];
  Rijm[1][1][2] = Rijmarray[6];
  Rijm[1][2][1] = Rijmarray[6];
  Rijm[2][1][1] = Rijmarray[6];
  Rijm[0][2][2] = Rijmarray[7];
  Rijm[2][0][2] = Rijmarray[7];
  Rijm[2][2][0] = Rijmarray[7];
  Rijm[1][2][2] = Rijmarray[8];
  Rijm[2][1][2] = Rijmarray[8];
  Rijm[2][2][1] = Rijmarray[8];
  Rijm[0][1][2] = Rijmarray[9];
  Rijm[1][2][0] = Rijmarray[9];
  Rijm[2][0][1] = Rijmarray[9];
  Rijm[2][1][0] = Rijmarray[9];
  Rijm[1][0][2] = Rijmarray[9];
  Rijm[0][2][1] = Rijmarray[9];

  Rmpp[0] = Rijm[0][0][0] + Rijm[1][1][0] + Rijm[2][2][0];
  Rmpp[1] = Rijm[0][0][1] + Rijm[1][1][1] + Rijm[2][2][1];
  Rmpp[2] = Rijm[0][0][2] + Rijm[1][1][2] + Rijm[2][2][2];

  bxdl[2] = burgX * dly - burgY * dlx;
  bxdl[0] = burgY * dlz - burgZ * dly;
  bxdl[1] = burgZ * dlx - burgX * dlz;

  bxRmpp[2] = burgX * Rmpp[1] - burgY * Rmpp[0];
  bxRmpp[0] = burgY * Rmpp[2] - burgZ * Rmpp[1];
  bxRmpp[1] = burgZ * Rmpp[0] - burgX * Rmpp[2];

  /* optimization: unroll */

  Sigma[0][0] = -1.0 * (bxRmpp[0] * dl[0]);
  Sigma[0][1] = Sigma[1][0] = -0.5 * (bxRmpp[0] * dl[1] + bxRmpp[1] * dl[0]);
  Sigma[0][2] = Sigma[2][0] = -0.5 * (bxRmpp[0] * dl[2] + bxRmpp[2] * dl[0]);
  Sigma[1][1] = -1.0 * (bxRmpp[1] * dl[1]);
  Sigma[1][2] = Sigma[2][1] = -0.5 * (bxRmpp[1] * dl[2] + bxRmpp[2] * dl[1]);
  Sigma[2][2] = -1.0 * (bxRmpp[2] * dl[2]);

  Rijm[0][0][0] -= Rmpp[0];
  Rijm[0][0][1] -= Rmpp[1];
  Rijm[0][0][2] -= Rmpp[2];

  Rijm[1][1][0] -= Rmpp[0];
  Rijm[1][1][1] -= Rmpp[1];
  Rijm[1][1][2] -= Rmpp[2];

  Rijm[2][2][0] -= Rmpp[0];
  Rijm[2][2][1] -= Rmpp[1];
  Rijm[2][2][2] -= Rmpp[2];

  Sigma[0][0] = MUover4pi *
                (Sigma[0][0] +
                 onepois * (Rijm[0][0][0] * bxdl[0] + Rijm[0][0][1] * bxdl[1] +
                            Rijm[0][0][2] * bxdl[2]));
  Sigma[0][1] = MUover4pi *
                (Sigma[0][1] +
                 onepois * (Rijm[0][1][0] * bxdl[0] + Rijm[0][1][1] * bxdl[1] +
                            Rijm[0][1][2] * bxdl[2]));
  Sigma[0][2] = MUover4pi *
                (Sigma[0][2] +
                 onepois * (Rijm[0][2][0] * bxdl[0] + Rijm[0][2][1] * bxdl[1] +
                            Rijm[0][2][2] * bxdl[2]));
  Sigma[1][1] = MUover4pi *
                (Sigma[1][1] +
                 onepois * (Rijm[1][1][0] * bxdl[0] + Rijm[1][1][1] * bxdl[1] +
                            Rijm[1][1][2] * bxdl[2]));
  Sigma[1][2] = MUover4pi *
                (Sigma[1][2] +
                 onepois * (Rijm[1][2][0] * bxdl[0] + Rijm[1][2][1] * bxdl[1] +
                            Rijm[1][2][2] * bxdl[2]));
  Sigma[2][2] = MUover4pi *
                (Sigma[2][2] +
                 onepois * (Rijm[2][2][0] * bxdl[0] + Rijm[2][2][1] * bxdl[1] +
                            Rijm[2][2][2] * bxdl[2]));

  Sigma[1][0] = Sigma[0][1];
  Sigma[2][0] = Sigma[0][2];
  Sigma[2][1] = Sigma[1][2];
}
