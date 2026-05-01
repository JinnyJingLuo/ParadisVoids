/***************************************************************************
 *
 *      Module:         Decomp.h
 *      Description:    This contains structures and prototypes needed
 *                      in modules that deal with the generic functions
 *                      for creating, accessing and manipulating the
 *                      domain decomposition.
 *
 ***************************************************************************/
#ifndef _Decomp_h
#define _Decomp_h

/*
 *      The dynamic load balancing is done based on per-process load
 *      data.  Below are definitions for the various types of load
 *      information that might be used.
 */
#define DLB_USE_WALLCLK_TIME 0
#define DLB_USE_FORCECALC_COUNT 1

/*
 *      When allocating an RB decomposition structure, in some cases
 *      we want to allocate a new decomposition and initialize it from
 *      scratch, in other cases, we want to duplicate an existing
 *      decomposition.  Define values to indicate to the allocate
 *      function which behaviour is desired.
 */
#define ALLOC_NEW_DECOMP 0
#define ALLOC_DUPLICATE_DECOMP 1

/*
 *      The RSDecomp_t structure defines arrays which contain the domain
 *      boundaries for the entire problem based on the original (recursive
 *      sectioning) version of the domain decomposition. There are only
 *      nXdoms partitions in the x direction. For each such partition
 *      there are nYdoms partitions, not necessarily at the same place
 *      for different x partitions. Thus a 2D array is needed.  Similarly,
 *      a 3D array  is needed for the final partitioning in the z direction.
 *
 *      See comments in RSDecomp.c for more details.
 */
typedef struct {
  real8 *domBoundX;
  real8 **domBoundY;
  real8 ***domBoundZ;
} RSDecomp_t;

/*
 *      Include prototypes for the set of generic function calls
 *      used to obtain information related to the domain decomposition.
 */
void BroadcastDecomp(Home_t *home, void *decomp);
void DLBfreeOld(Home_t *home);
int FindCoordDomain(Home_t *home, const double &x, const double &y,
                    const double &z);
int FindCoordDomain(Home_t *home, Point *poPoint);
void FreeDecomp(Home_t *home, void *decomp);
void GetAllDecompBounds(Home_t *home, real8 **decompBounds, int *numValues);
void GetCellDomainList(Home_t *home, int cellID, int *domCount, int **domList);
void GetLocalDomainBounds(Home_t *home, void *decomp);
void Rebalance(Home_t *home);
void UniformDecomp(Home_t *home, void **decomp);
RSDecomp_t *PointBasedDecomp(Home_t *home, InData_t *poInputData);
RSDecomp_t *PointBasedStressDecomp(Home_t *home, InData_t *poInputData);
void WriteDecompBounds(Home_t *home, FILE *fp);

#endif
