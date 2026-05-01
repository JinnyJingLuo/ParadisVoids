#include "Home.h"
#include "Util.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <stdio.h>
#include <math.h>

/*---------------------------------------------------------------------------
 *
 *      Function:    Tecplot
 *      Description: Write segment data in tecplot format... Masato
 *
 *      Args:
 *          baseFileName     Base name of the plot file.  Plot data
 *                           will be written to 1 or more file segments
 *                           named <baseFileName>.n
 *          ioGroup          I/O group number associated with this domain
 *          firstInGroup     1 if this domain is the first processor in
 *                           its I/O group, zero otherwise.
 *          writePrologue    1 if this process should write all needed
 *                           headers and do any initialization associated
 *                           with the plot file, zero otherwise
 *          writeEpilogue    1 if this process should write all needed
 *                           trailers and do any terminal processing
 *                           associated with the plot file, zero otherwise
 *          numSegs          count of all segments in the problem space
 *
 *-------------------------------------------------------------------------*/
void Tecplot(Home_t *home, char *baseFileName, int ioGroup, int firstInGroup,
             int writePrologue, int writeEpilogue, int numSegs) {
  int i, j, thisDomain;
  int newNodeKeyPtr;
  int btype;
  real8 x, y, z;
  real8 x2, y2, z2;
  real8 x2o, y2o, z2o;
  real8 prx, pry, prz;
  real8 vx, vy, vz;
  real8 bX, bY, bZ;
  real8 Lx, Ly, Lz;
  char fileName[256];
  Node_t *node, *nbrNode;
  Param_t *param;
  FILE *fp;
  struct stat statbuf;

  param = home->param;
  thisDomain = home->myDomain;

  Lx = param->Lx;
  Ly = param->Ly;
  Lz = param->Lz;

  /*
   *      Set data file name.  Only append a sequence number to
   *      the data file name if the data is to be spread across
   *      multiple files.
   */
  if (param->numIOGroups == 1) {
    snprintf(fileName, sizeof(fileName), "%s/%s.dat", DIR_TECPLOT,
             baseFileName);
  } else {
    snprintf(fileName, sizeof(fileName), "%s/%s.%d", DIR_TECPLOT, baseFileName,
             ioGroup);
  }

#ifdef PARALLEL
#ifdef DO_IO_TO_NFS
  /*
   *      It appears that when multiple processes on different hosts
   *      write to the same NFS-mounted file (even if access is
   *      sequentialized), consistency problems can arise resulting
   *      in corrupted data files.  Explicitly doing a 'stat' of the
   *      file on each process immediately before opening it *seems*
   *      to clear up the problem.
   */
  memset(&statbuf, 0, sizeof(statbuf));
  (void)stat(fileName, &statbuf);
#endif
#endif

  /*
   *      If this process is the first member of its I/O group
   *      it needs to create the output file and do any intializations
   *      needed.
   */
  if (firstInGroup) {
    /*
     *          First task in the I/O group must open the data file for writing
     *          to overwrite any existing file of the same name.
     */
    if ((fp = fopen(fileName, "w")) == (FILE *)NULL) {
      Fatal("tec_plot: Open error %d on %s\n", errno, fileName);
    }

    if (writePrologue) {
      printf(" +++ Writing Tecplot file(s) %s\n", baseFileName);

      fprintf(fp, "variables = X,Y,Z,V1,V2,V3,V4,V5,V6,V7\n");
      fprintf(fp, "zone i = %d  F=POINT\n", 2 * numSegs);
    }
  } else {
    /*
     *          Any task NOT first in its I/O group must open the data file
     *          in an append mode so everything gets added to the end of
     *          the file.
     */
    if ((fp = fopen(fileName, "a")) == (FILE *)NULL) {
      Fatal("tec_plot: Open error %d on %s\n", errno, fileName);
    }
  }

  /*
   *      Generate the plot data for the segments associated with this
   *      domain's data
   */
  newNodeKeyPtr = home->newNodeKeyPtr;
  double dTolerance = 1.0E-6;
  for (i = 0; i < newNodeKeyPtr; i++) {
    node = home->nodeKeys[i];
    if (node == NULL)
      continue;
    x = node->x;
    y = node->y;
    z = node->z;

    for (j = 0; j < node->numNbrs; j++) {
      // do not add a segment twice, use the positive nonzero Burgers vector
      // component criterion
      bX = node->burgX[j];
      bY = node->burgY[j];
      bZ = node->burgZ[j];
      if (fabs(bX) < dTolerance) {
        if (fabs(bY) < dTolerance) {
          if (bZ <= dTolerance) {
            continue;
          }
        } else if (bY < 0.0) {
          continue;
        }
      } else if (bX < 0.0) {
        continue;
      }

      if (fabs(bX * bY * bZ) < dTolerance) {
        if (fabs(bX * bY) > 0.0 || fabs(bY * bZ) > 0.0 || fabs(bZ * bX) > 0.0) {
          // a [110] type Burgers vector
          btype = 10;
        } else {
          // a [100] type Burgers vector
          btype = 0;
        }
      } else {
        if (bX * bY * bZ < 0.0) {
          bX *= -1;
          bY *= -1;
          bZ *= -1;
        }
        if (fabs(fabs(bX) - fabs(bY)) > dTolerance ||
            fabs(fabs(bY) - fabs(bZ)) > dTolerance ||
            fabs(fabs(bZ) - fabs(bX)) > dTolerance) {
          // a non [111],[110],[100] or [000] type vector
          btype = 20;
        } else {
          // a [111] type Burgers vector
          if (bY < 0 && bZ < 0)
            btype = 1;
          if (bZ < 0 && bX < 0)
            btype = 2;
          if (bX < 0 && bY < 0)
            btype = 3;
          if (bX > 0 && bY > 0 && bZ > 0)
            btype = 4;
        }
      }

      nbrNode = GetNeighborNode(home, node, j);

      if (nbrNode == NULL) {
        printf("WARNING: Neighbor not found at %s line %d\n", __FILE__,
               __LINE__);
        continue;
      }

      x2o = nbrNode->x;
      y2o = nbrNode->y;
      z2o = nbrNode->z;

      x2 = x2o;
      y2 = y2o;
      z2 = z2o;
      GetNearestImage(home->param, x, y, z, &x2, &y2, &z2);
      prx = x2 - x;
      pry = y2 - y;
      prz = z2 - z;
      fprintf(
          fp, "%15.10f %15.10f %15.10f %15.10f %15.10f %15.10f %15.10f %15.10f "
              "%15.10f %d\n",
          x, y, z, prx, pry, prz, node->nx[j], node->ny[j], node->nz[j], btype);

      x2 = x;
      y2 = y;
      z2 = z;
      GetNearestImage(home->param, x2o, y2o, z2o, &x2, &y2, &z2);
      prx = x2 - x2o;
      pry = y2 - y2o;
      prz = z2 - z2o;
      fprintf(fp,
              "%15.10f %15.10f %15.10f %15.10f %15.10f %15.10f %15.10f %15.10f "
              "%15.10f %d\n",
              x2o, y2o, z2o, -prx, -pry, -prz, node->nx[j], node->ny[j],
              node->nz[j], btype);
    }
  }

  /*
   *      Handle anything that needs to be done after all data has been
   *      written to the file
   */
  if (writeEpilogue) {
    fprintf(fp, "\n");
  }

  fclose(fp);
}
