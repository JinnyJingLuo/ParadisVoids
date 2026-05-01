#ifndef _FORCE_H
#define _FORCE_H
/***************************************************************************
 *
 *	Module:		Force.h
 *	Description:	This header is primarily for prototypes of various
 *                      functions used in calculating or updating
 *                      forces or stresses.
 *
 ***************************************************************************/
#include "Home.h"
#include "ParadisExternalLoadServer.h"

void AddtoArmForce(Node_t *node, int arm, real8 f[3]);
void AddtoNodeForce(Node_t *node, real8 f[3]);
void ComputeForces(Home_t *home, Node_t *node1, Node_t *node2, Node_t *node3,
                   Node_t *node4, real8 *f1, real8 *f2, real8 *f3, real8 *f4);
void ComputeSegSigbRem(Home_t *home, int reqType);
void deWitInteraction(real8 MU, real8 NU, real8 Sigma[][3], real8 px, real8 py,
                      real8 pz, real8 tp1, real8 tp2, real8 tp3, real8 burgX,
                      real8 burgY, real8 burgZ, real8 *cntx, real8 *cnty,
                      real8 *cntz);
void dSegImgStress(Home_t *home, real8 Sigma[][3], real8 px, real8 py, real8 pz,
                   real8 dlx, real8 dly, real8 dlz, real8 burgX, real8 burgY,
                   real8 burgZ, real8 rx, real8 ry, real8 rz, int pbc);
void GetFieldPointStress(Home_t *home, real8 x, real8 y, real8 z,
                         real8 totStress[3][3]);
void LocalSegForces(Home_t *home,
                    ParadisExternalLoadServer *poExternalLoadServer,
                    ParadisPrecipitateServer *poPrecipitateServer);
void NodeForce(Home_t *home, ParadisExternalLoadServer *poExternalLoadServer,
               ParadisPrecipitateServer *poPrecipitateServer);
void OsmoticForce(Home_t *home, real8 x1, real8 y1, real8 z1, real8 x2,
                  real8 y2, real8 z2, real8 bx, real8 by, real8 bz, real8 f1[3],
                  real8 f2[3]);
void PKForce(real8 sigb[3], real8 X1, real8 Y1, real8 Z1, real8 X2, real8 Y2,
             real8 Z2, real8 f1[3], real8 f2[3]);
void ReevaluateForces(Home_t *home);

#ifdef _FEM
void SegmentStress(real8 MU, real8 NU, real8 burgX, real8 burgY, real8 burgZ,
                   real8 xA, real8 yA, real8 zA, real8 xB, real8 yB, real8 zB,
                   real8 x0, real8 y0, real8 z0, real8 a, real8 Sigma[3][3]);
#endif

void SegSegForceIsotropic(real8 p1x, real8 p1y, real8 p1z, real8 p2x, real8 p2y,
                          real8 p2z, real8 p3x, real8 p3y, real8 p3z, real8 p4x,
                          real8 p4y, real8 p4z, real8 bpx, real8 bpy, real8 bpz,
                          real8 bx, real8 by, real8 bz, real8 a, real8 MU,
                          real8 NU, int seg12Local, int seg34Local, real8 *fp1x,
                          real8 *fp1y, real8 *fp1z, real8 *fp2x, real8 *fp2y,
                          real8 *fp2z, real8 *fp3x, real8 *fp3y, real8 *fp3z,
                          real8 *fp4x, real8 *fp4y, real8 *fp4z);
void SegSegForce(real8 p1x, real8 p1y, real8 p1z, real8 p2x, real8 p2y,
                 real8 p2z, real8 p3x, real8 p3y, real8 p3z, real8 p4x,
                 real8 p4y, real8 p4z, real8 bpx, real8 bpy, real8 bpz,
                 real8 bx, real8 by, real8 bz, real8 a, real8 MU, real8 NU,
                 int seg12Local, int seg34Local, real8 *fp1x, real8 *fp1y,
                 real8 *fp1z, real8 *fp2x, real8 *fp2y, real8 *fp2z,
                 real8 *fp3x, real8 *fp3y, real8 *fp3z, real8 *fp4x,
                 real8 *fp4y, real8 *fp4z);
void SelfForce(int coreOnly, real8 MU, real8 NU, real8 bx, real8 by, real8 bz,
               real8 x1, real8 y1, real8 z1, real8 x2, real8 y2, real8 z2,
               real8 a, real8 Ecore, real8 f1[3], real8 f2[3]);
void SemiInfiniteSegSegForce(real8 p1x, real8 p1y, real8 p1z, real8 p2x,
                             real8 p2y, real8 p2z, real8 p3x, real8 p3y,
                             real8 p3z, real8 p4x, real8 p4y, real8 p4z,
                             real8 bpx, real8 bpy, real8 bpz, real8 bx,
                             real8 by, real8 bz, real8 a, real8 MU, real8 NU,
                             real8 *fp1x, real8 *fp1y, real8 *fp1z, real8 *fp3x,
                             real8 *fp3y, real8 *fp3z, real8 *fp4x, real8 *fp4y,
                             real8 *fp4z);
void ZeroNodeForces(Home_t *home, int reqType);

#endif /* ifndef _Force_h */
