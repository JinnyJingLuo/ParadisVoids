#include "Home.h"
#include "TwinPlaneCrossSlip.h"
#include "ParadisBlockSurface.h"
#include <cassert>
using namespace std;

void TwinPlaneCrossSlip(Home_t *home) {
  if (home->param->EnableTwinPlaneCrossSlip == 0) {
    return;
  }
  // qjiao: interface cross-slip
  // 2019/09/12: updated if..else algorithm for determine the node->csState (v3)
  // 2019/09/14: handle dislocations not able to cross-slip (v4)
  Node_t *poNbr1, *poNbr2;
  Node_t *node;
  double dTol = 1E-6; // set a tolerance for geometrical detection
  double len;
  double A, B, C, D;
  double X_1, Y_1, Z_1;
  double daTemp[3] = {0.0, 0.0, 0.0};

  // remember: you also need to change the definition of cs plane 1.in
  // crossslipserver 2.in s-s sgmentcollision 3.in hingejointcollision
  A = home->param->A;
  B = home->param->B;
  C = home->param->C;
  X_1 = home->param->X_1;
  Y_1 = home->param->Y_1;
  Z_1 = home->param->Z_1;
  D = -(A * X_1 + B * Y_1 + C * Z_1);
  // fprintf(stderr, "Warning\n");
  // printf("Warning!\n");
  len = sqrt(pow(A, 2.0) + pow(B, 2.0) + pow(C, 2.0));
  real8 segPlaneX = A / len;
  real8 segPlaneY = B / len;
  real8 segPlaneZ = C / len;
  real8 dx, dy, dz;
  real8 lenDeduc;
  // fprintf(stderr, "Length = %e\n", len);
  Point pOld, pNew, oIntersectionPoint;
  Point pOld_ex; // junjie: extrapolation
  Line nodeTrace;
  Vector poPlaneNorm, crossSlipPlaneNorm;
  Plane poPlane, crossSlipPlane;
  Vector disp;
  Vector oNormal1, oNormal2, oline;
  bool IsIntersection;

  crossSlipPlane.Set(Vector(A, B, C), Point(X_1, Y_1, Z_1));

  for (int i = 0; i < home->newNodeKeyPtr; i++) {
    if ((node = home->nodeKeys[i]) == NULL)
      continue;
    if (node->constraint == PINNED_NODE)
      continue;
    if (node->myTag.domainID != home->myDomain)
      continue;
    node->csState = ON_SLIP_PLANE;
    pOld.Set(node->oldx, node->oldy, node->oldz);
    pNew.Set(node->x, node->y, node->z);
    // case that node not glide on the cross-slip plane
    if ((crossSlipPlane.GetPointDistance(pOld) > dTol) ||
        (crossSlipPlane.GetPointDistance(pNew) > dTol)) {
      if (crossSlipPlane.GetSegmentIntersection(pOld, pNew,
                                                oIntersectionPoint)) {
        //	fprintf(stderr, "Xing = 1, ");
        // cross-slip determination
        IsIntersection = true;
        for (int j = 0; j < node->numNbrs; ++j) {
          if (abs(A * node->burgX[j] + B * node->burgY[j] +
                  C * node->burgZ[j]) > dTol)
            IsIntersection = false;
        }
        if (IsIntersection) {
          //		fprintf(stderr, "Op = RP, ");
          daTemp[0] = oIntersectionPoint.GetX();
          daTemp[1] = oIntersectionPoint.GetY();
          daTemp[2] = oIntersectionPoint.GetZ();
          RepositionNode(home, daTemp, &(node->myTag), 1);
          node->csState = AT_INTERSECTION;
          // continue;  //need this continue becase next if case if
          // AT_INTERSECTION
        } else {
          double daPos[3] = {0.0, 0.0, 0.0};
          daPos[0] = node->oldx;
          daPos[1] = node->oldy;
          daPos[2] = node->oldz;
          RepositionNode(home, daPos, &(node->myTag), 1);
        }
      } else if (crossSlipPlane.GetPointDistance(pOld) <= dTol &&
                 crossSlipPlane.GetPointDistance(pNew) > dTol) {
        dx = node->x - node->oldx;
        dy = node->y - node->oldy;
        dz = node->z - node->oldz;
        disp = Vector(dx, dy, dz);
        oNormal1 = Vector(segPlaneX, segPlaneY, segPlaneZ);
        oNormal2 = Vector(node->nx[0], node->ny[0], node->nz[0]);
        oline = oNormal1 ^ oNormal2;
        lenDeduc = disp * oline;
        disp = oline * lenDeduc;
        daTemp[0] = node->oldx + disp.GetX();
        daTemp[1] = node->oldy + disp.GetY();
        daTemp[2] = node->oldz + disp.GetZ();
        RepositionNode(home, daTemp, &(node->myTag), 1);
      } else if (crossSlipPlane.GetPointDistance(pOld) > dTol &&
                 crossSlipPlane.GetPointDistance(pNew) <= dTol) {
        double daPos[3] = {0.0, 0.0, 0.0};
        daPos[0] = node->oldx;
        daPos[1] = node->oldy;
        daPos[2] = node->oldz;
        RepositionNode(home, daPos, &(node->myTag), 1);
        printf("dangerous!\n");
      } else if (crossSlipPlane.GetPointDistance(pOld) > dTol &&
                 crossSlipPlane.GetPointDistance(pNew) > dTol)
        ;
      else
        assert(false);
    } else {
      Vector b1;
      bool bHave112, bHave110;
      bHave112 = false;
      bHave110 = false;
      for (int j = 0; j < node->numNbrs; ++j) {
        b1 = Vector(node->burgX[j], node->burgY[j], node->burgZ[j]);
        if (ParadisCrossSlipServer::Is112Vector(b1)) {
          bHave112 = true;
          break;
        }
      }
      for (int j = 0; j < node->numNbrs; ++j) {
        b1 = Vector(node->burgX[j], node->burgY[j], node->burgZ[j]);
        if (ParadisCrossSlipServer::Is110Vector(b1)) {
          bHave110 = true;
          break;
        }
      }
      if (bHave110 && bHave112) {
        node->csState = AT_INTERSECTION;
        continue;
      } else if (!bHave110 && bHave112) {
        node->csState = ON_CROSS_SLIP_PLANE;
        continue;
      } else if (bHave110 && !bHave112) {
        node->csState = AT_INTERSECTION;
        continue;
      }
      assert(false);
    }
  }
  /*
  for (int i = 0; i < home->newNodeKeyPtr; i++) {
      if ((node = home->nodeKeys[i]) == NULL) continue;
      if (node->constraint == PINNED_NODE) continue;
      if (node->myTag.domainID != home->myDomain) continue;
      if (node->csState == AT_INTERSECTION) {
          if (node->constraint >=SURFACE_NODE) {
              pNew.Set(node->x, node->y, node->z);
              if (crossSlipPlane.GetPointDistance(pNew) > dTol) {
                  daTemp[0] = node->oldx;
                  daTemp[1] = node->oldy;
                  daTemp[2] = node->oldz;
                  RepositionNode(home, daTemp, &(node->myTag), 1);
              }
          }
          else
          {
              unsigned int iConstraintType = 0;
              Vector oConstraintVector;
              ParadisSurface::GetNodeDynamicConstraint(node,iConstraintType,oConstraintVector);
              if(iConstraintType >= 2) continue;
              else {
                  dx = node->x - node->oldx;
                  dy = node->y - node->oldy;
                  dz = node->z - node->oldz;
                  disp = Vector(dx, dy, dz);
                  oNormal1 = Vector(segPlaneX, segPlaneY, segPlaneZ);
                  oNormal2 = Vector(node->nx[0], node->ny[0], node->nz[0]);
                  oline = oNormal1 ^ oNormal2;
                  lenDeduc = disp * oline;
                  disp = oline * lenDeduc;
                  pNew.Set(node->x, node->y, node->z);

                  if (crossSlipPlane.GetPointDistance(pNew) > dTol) {
                      daTemp[0] = node->oldx + disp.GetX();
                      daTemp[1] = node->oldy + disp.GetY();
                      daTemp[2] = node->oldz + disp.GetZ();
                      RepositionNode(home, daTemp, &(node->myTag), 1);
                  }
              }
          }
      }
  }
  */
  for (int i = 0; i < home->newNodeKeyPtr; i++) {
    if ((node = home->nodeKeys[i]) == NULL)
      continue;
    if (node->constraint == PINNED_NODE)
      continue;
    if (node->myTag.domainID != home->myDomain)
      continue;
    if (node->csState == AT_INTERSECTION) {
      if (node->constraint >=SURFACE_NODE && node->numNbrs == 1) {
        poNbr1 = GetNeighborNode(home, node, 0);
        if (poNbr1 == (Node_t *)NULL) {
          Fatal("Neighbor not found at %s line %d\n", __FILE__, __LINE__);
        }
        if (poNbr1->numNbrs == 3) {
          SurfaceNodeDissociation(home, node);
        }
      }
      if (node->numNbrs != 2)
        continue;
      poNbr1 = GetNeighborNode(home, node, 0);
      poNbr2 = GetNeighborNode(home, node, 1);
      if (poNbr1 == (Node_t *)NULL || poNbr2 == (Node_t *)NULL) {
        Fatal("Neighbor not found at %s line %d\n", __FILE__, __LINE__);
      }
      if (poNbr1->myTag.domainID != home->myDomain ||
          poNbr2->myTag.domainID != home->myDomain)
        continue;

      /*            if (poNbr1->myTag.domainID != home->myDomain) {
                      pOld.Set(poNbr1->oldx, poNbr1->oldy, poNbr1->oldz);
                      pNew.Set(poNbr1->x, poNbr1->y, poNbr1->z);
                      // it could cross the interface if it was already
         repositioned at the interface
                      // from previous step
                      if (crossSlipPlane.GetSegmentIntersection(pOld, pNew,
         oIntersectionPoint) || (crossSlipPlane.GetPointDistance(pOld) < dTol))
         { poNbr1->csState = AT_INTERSECTION;
                      }
                              }
                  if (poNbr2->myTag.domainID != home->myDomain) {
                      pOld.Set(poNbr2->oldx, poNbr2->oldy, poNbr2->oldz);
                      pNew.Set(poNbr2->x, poNbr2->y, poNbr2->z);
                      // it could cross the interface if it was already
         repositioned at the interface
                      // from previous step
                      if (crossSlipPlane.GetSegmentIntersection(pOld, pNew,
         oIntersectionPoint) || (crossSlipPlane.GetPointDistance(pOld) < dTol))
         { poNbr2->csState = AT_INTERSECTION;
                      }
                              }*/

      if (poNbr1->csState != AT_INTERSECTION ||
          poNbr2->csState != AT_INTERSECTION)
        continue;
      Node_t *poend1, *poend2, *pocenter;
      poend1 = poNbr1;
      poend2 = poNbr2;
      pocenter = node;
      bool bEnd1found = false;
      bool bEnd2found = false;
      int numNodeIntersection = 3;
      // see if poend1 is the end of the intersection line
      while (bEnd1found == false) {
        if (poend1->numNbrs == 2) {
          poNbr1 = GetNeighborNode(home, poend1, 0);
          poNbr2 = GetNeighborNode(home, poend1, 1);
          if (poNbr1 == (Node_t *)NULL || poNbr2 == (Node_t *)NULL) {
            Fatal("Neighbor not found at %s line %d\n", __FILE__, __LINE__);
          }
          if (poNbr1->myTag.domainID != home->myDomain ||
              poNbr2->myTag.domainID != home->myDomain) {
            bEnd1found = true;
            continue;
          }

          if (poNbr1->csState == AT_INTERSECTION &&
              poNbr2->csState == AT_INTERSECTION) {
            numNodeIntersection++;
            if (poNbr1 != pocenter) {
              pocenter = poend1;
              poend1 = poNbr1;
            } else {
              pocenter = poend1;
              poend1 = poNbr2;
            }
          } else
            bEnd1found = true;
        } else
          bEnd1found = true;
      }
      // see if poend2 is the end of the intersection line
      pocenter = node;
      while (bEnd2found == false) {
        if (poend2->numNbrs == 2) {
          poNbr1 = GetNeighborNode(home, poend2, 0);
          poNbr2 = GetNeighborNode(home, poend2, 1);
          if (poNbr1 == (Node_t *)NULL || poNbr2 == (Node_t *)NULL) {
            Fatal("Neighbor not found at %s line %d\n", __FILE__, __LINE__);
          }
          if (poNbr1->myTag.domainID != home->myDomain ||
              poNbr2->myTag.domainID != home->myDomain) {
            bEnd2found = true;
            continue;
          }

          if (poNbr1->csState == AT_INTERSECTION &&
              poNbr2->csState == AT_INTERSECTION) {
            numNodeIntersection++;
            if (poNbr1 != pocenter) {
              pocenter = poend2;
              poend2 = poNbr1;
            } else {
              pocenter = poend2;
              poend2 = poNbr2;
            }
          } else
            bEnd2found = true;
        } else
          bEnd2found = true;
      }
      Node_t *newNode = NULL;
      Dissociate(home, pocenter, numNodeIntersection, poend1, poend2, 1,
                 newNode);
    }
  }
}

/*---------------------------------------------------------------------------
 *
 *      Function:       DissociateFourArmsNode
 *      Description:    Merge two 3-arm nodes into a four-arm node when they
 *                      get close enough to each other on the intersections
 *
 *      junjie
 *
 *-------------------------------------------------------------------------*/
void MakeFourArmsNode(Home_t *home) {
  if (home->param->EnableTwinPlaneCrossSlip == 0) {
    return;
  }
  Node_t *node;
  for (int k = 0; k < home->newNodeKeyPtr; k++) {
    if ((node = home->nodeKeys[k]) == NULL)
      continue;
    if (node->constraint == PINNED_NODE)
      continue;
    if (node->myTag.domainID != home->myDomain)
      continue;
    if (node->csState != AT_INTERSECTION)
      continue;
    if (node->numNbrs != 3)
      continue;
    // this node must have 2 <112> vectors and 1 <110> vectors, find them
    Vector b1, b2, b3;
    real8 bx1, by1, bz1;
    real8 bx2, by2, bz2;
    real8 bx3, by3, bz3;
    real8 nx, ny, nz;
    real8 dx, dy, dz;
    Node_t *poNbr1, *poNbr2, *poNbr3, *poNbr4, *poNbr5;
    int found110;
    int found112;
    real8 dTol = 1E-6;
    double A, B, C, len;
    A = home->param->A;
    B = home->param->B;
    C = home->param->C;
    len = sqrt(pow(A, 2.0) + pow(B, 2.0) + pow(C, 2.0));
    real8 segPlaneX = A / len;
    real8 segPlaneY = B / len;
    real8 segPlaneZ = C / len;
    found110 = 0;
    found112 = 0;
    for (int i = 0; i < node->numNbrs; ++i) {
      poNbr1 = GetNeighborNode(home, node, i);
      if (poNbr1 == (Node_t *)NULL) {
        Fatal("Neighbor not found at %s line %d\n", __FILE__, __LINE__);
      }
      GetBurgersVectorNormal(home, node, poNbr1, &bx1, &by1, &bz1, &nx, &ny,
                             &nz);
      b1 = Vector(bx1, by1, bz1);
      if (ParadisCrossSlipServer::Is110Vector(b1)) {
        found110 = 1;
        if (pow(node->nx[i] - segPlaneX, 2) + pow(node->ny[i] - segPlaneY, 2) +
                pow(node->nz[i] - segPlaneZ, 2) <
            dTol)
          continue;
        if (found112 == 0) {
          poNbr2 = GetNeighborNode(home, node, 1);
          poNbr3 = GetNeighborNode(home, node, 2);
          if (poNbr2 == (Node_t *)NULL || poNbr3 == (Node_t *)NULL) {
            Fatal("Neighbor not found at %s line %d\n", __FILE__, __LINE__);
          }
        } else if (found112 == 1) {
          poNbr3 = GetNeighborNode(home, node, 2);
          if (poNbr3 == (Node_t *)NULL) {
            Fatal("Neighbor not found at %s line %d\n", __FILE__, __LINE__);
          }
        }
        break;
      } else if (found112 == 0) {
        poNbr2 = GetNeighborNode(home, node, i);
        if (poNbr2 == (Node_t *)NULL) {
          Fatal("Neighbor not found at %s line %d\n", __FILE__, __LINE__);
        }
        ++found112;
      } else if (found112 == 1) {
        poNbr3 = GetNeighborNode(home, node, i);
        if (poNbr3 == (Node_t *)NULL) {
          Fatal("Neighbor not found at %s line %d\n", __FILE__, __LINE__);
        }
        ++found112;
      }
    }
    if (found110 != 1)
      continue;
    GetBurgersVectorNormal(home, node, poNbr2, &bx2, &by2, &bz2, &nx, &ny, &nz);
    b2 = Vector(bx2, by2, bz2);
    GetBurgersVectorNormal(home, node, poNbr3, &bx3, &by3, &bz3, &nx, &ny, &nz);
    b3 = Vector(bx3, by3, bz3);
    if (!ParadisCrossSlipServer::Is112Vector(b2) ||
        !ParadisCrossSlipServer::Is112Vector(b3))
      continue;

    if (poNbr1->myTag.domainID != home->myDomain)
      continue;
    if (poNbr1->csState != AT_INTERSECTION)
      continue;
    if (poNbr1->numNbrs != 3)
      continue;

    dx = node->x - poNbr1->x;
    dy = node->y - poNbr1->y;
    dz = node->z - poNbr1->z;
    real8 dNewPosition[3];
    Node_t *poTargetNode = NULL;
    dNewPosition[0] = (node->x + poNbr1->x) / 2.0;
    dNewPosition[1] = (node->y + poNbr1->y) / 2.0;
    dNewPosition[2] = (node->z + poNbr1->z) / 2.0;
    if (sqrt(pow(dx, 2) + pow(dy, 2) + pow(dz, 2)) <= home->param->minSeg) {
      // Merge node and poNbr1
      MergeNode(home, OPCLASS_COLLISION, node, poNbr1, dNewPosition,
                &poTargetNode, 1);
    }
  }
}

/*---------------------------------------------------------------------------
 *
 *      Function:       DissociateFourArmsNode
 *      Description:    Dissociate the four arms node on the intersections
 *
 *      junjie
 *
 *-------------------------------------------------------------------------*/
void DissociateFourArmsNode(Home_t *home) {
  if (home->param->EnableTwinPlaneCrossSlip == 0) {
    return;
  }
  MakeFourArmsNode(home);
  real8 dTol = 1E-6;
  Node_t *node;
  Node_t *poNbr1, *poNbr2, *poNbr3, *poNbr4;
  real8 bx1, by1, bz1;
  real8 bx2, by2, bz2;
  real8 bx3, by3, bz3;
  real8 bx4, by4, bz4;
  Vector b1, b2, b3, b4;
  real8 bx, by, bz;
  real8 nx, ny, nz;
  real8 segPlaneX, segPlaneY, segPlaneZ;
  real8 A, B, C, len;
  bool ThisNodeShouldBeDissociated = false;
  A = home->param->A;
  B = home->param->B;
  C = home->param->C;
  len = sqrt(pow(A, 2.0) + pow(B, 2.0) + pow(C, 2.0));
  segPlaneX = A / len;
  segPlaneY = B / len;
  segPlaneZ = C / len;

  for (int i = 0; i < home->newNodeKeyPtr; i++) {
    if ((node = home->nodeKeys[i]) == NULL)
      continue;
    if (node->constraint == PINNED_NODE)
      continue;
    if (node->myTag.domainID != home->myDomain)
      continue;
    if (node->csState == ON_SLIP_PLANE)
      continue;
    // if (node->csState != AT_INTERSECTION) continue;
    if (node->numNbrs != 4)
      continue;
    if (node->flags & NO_COLLISIONS)
      continue;
    ThisNodeShouldBeDissociated = true;
    poNbr1 = GetNeighborNode(home, node, 0);
    poNbr2 = GetNeighborNode(home, node, 1);
    poNbr3 = GetNeighborNode(home, node, 2);
    poNbr4 = GetNeighborNode(home, node, 3);
    GetBurgersVectorNormal(home, node, poNbr1, &bx1, &by1, &bz1, &nx, &ny, &nz);
    b1 = Vector(bx1, by1, bz1);
    GetBurgersVectorNormal(home, node, poNbr2, &bx2, &by2, &bz2, &nx, &ny, &nz);
    b2 = Vector(bx2, by2, bz2);
    GetBurgersVectorNormal(home, node, poNbr3, &bx3, &by3, &bz3, &nx, &ny, &nz);
    b3 = Vector(bx3, by3, bz3);
    GetBurgersVectorNormal(home, node, poNbr4, &bx4, &by4, &bz4, &nx, &ny, &nz);
    b4 = Vector(bx4, by4, bz4);
    if (!ParadisCrossSlipServer::Is112Vector(b1) ||
        !ParadisCrossSlipServer::Is112Vector(b2) ||
        !ParadisCrossSlipServer::Is112Vector(b3) ||
        !ParadisCrossSlipServer::Is112Vector(b4)) {
      continue;
    }
    // now need to find out which 2 belong to the same side and the other 2
    // belong to the other side
    poNbr1 = GetNeighborNode(home, node, 0);
    if (poNbr1 == (Node_t *)NULL) {
      Fatal("Neighbor not found at %s line %d\n", __FILE__, __LINE__);
    }
    GetBurgersVectorNormal(home, node, poNbr1, &bx1, &by1, &bz1, &nx, &ny, &nz);
    for (int j = 1; j < node->numNbrs; ++j) {
      poNbr2 = GetNeighborNode(home, node, j);
      if (poNbr2 == (Node_t *)NULL) {
        Fatal("Neighbor not found at %s line %d\n", __FILE__, __LINE__);
      }
      GetBurgersVectorNormal(home, node, poNbr2, &bx, &by, &bz, &nx, &ny, &nz);
      if (pow(bx + bx1, 2) + pow(by + by1, 2) + pow(bz + bz1, 2) < dTol) {
        for (int k = 1; k < node->numNbrs; ++k) {
          if (k != j) {
            poNbr3 = GetNeighborNode(home, node, k);
            if (poNbr3 == (Node_t *)NULL) {
              Fatal("Neighbor not found at %s line %d\n", __FILE__, __LINE__);
            }
            GetBurgersVectorNormal(home, node, poNbr3, &bx3, &by3, &bz3, &nx,
                                   &ny, &nz);
            for (int l = 1; l < node->numNbrs; ++l) {
              if (l != j && l != k) {
                poNbr4 = GetNeighborNode(home, node, l);
                if (poNbr4 == (Node_t *)NULL) {
                  Fatal("Neighbor not found at %s line %d\n", __FILE__,
                        __LINE__);
                }
                GetBurgersVectorNormal(home, node, poNbr4, &bx4, &by4, &bz4,
                                       &nx, &ny, &nz);
                if (pow(bx3 + bx4, 2) + pow(by3 + by4, 2) + pow(bz3 + bz4, 2) >
                    dTol)
                  ThisNodeShouldBeDissociated = false;
                break;
              }
            }
            break;
          }
        }
        break;
      } else if (j == node->numNbrs - 1)
        ThisNodeShouldBeDissociated = false;
    }
    if (ThisNodeShouldBeDissociated == false)
      continue;
    Node_t *newNode;
    // now 1 & 2 are on the same side, 3 & 4 are on the same side
    newNode = GetNewNativeNode(home);
    FreeNodeArms(newNode);
    newNode->constraint = UNCONSTRAINED;
    newNode->dNx = node->dNx;
    newNode->dNy = node->dNy;
    newNode->dNz = node->dNz;
    newNode->oldx = node->x;
    newNode->oldy = node->y;
    newNode->oldz = node->z;
    newNode->x = (poNbr1->x + poNbr2->x) / 2;
    newNode->y = (poNbr1->y + poNbr2->y) / 2;
    newNode->z = (poNbr1->z + poNbr2->z) / 2;
    newNode->csState = ON_CROSS_SLIP_PLANE;

    if (1) {
      AddOp(home, SPLIT_NODE, node->myTag.domainID, node->myTag.index,
            newNode->myTag.domainID, newNode->myTag.index, -1,
            -1,            /* 3rd node tag unneeded */
            0.0, 0.0, 0.0, /* burgers vector unneeded*/
            newNode->x, newNode->y, newNode->z, 0.0, 0.0,
            0.0); // normal not needed
    }

    node->x = (poNbr3->x + poNbr4->x) / 2;
    node->y = (poNbr3->y + poNbr4->y) / 2;
    node->z = (poNbr3->z + poNbr4->z) / 2;
    node->csState = ON_CROSS_SLIP_PLANE;

    if (1) {
      AddOp(home, RESET_COORD, node->myTag.domainID, node->myTag.index, -1,
            -1,            /* 2nd tag unneeded */
            -1, -1,        /* 3rd node tag unneeded */
            0.0, 0.0, 0.0, /* burgers vector unneeded */
            node->x, node->y, node->z, 0.0, 0.0,
            0.0); /* plane normal not needed */
    }

    // connect
    int iChainID = 0;
    // junjie: Can I insert an arm to a ghost node???
    InsertArm(home, newNode, &poNbr1->myTag, bx1, by1, bz1, segPlaneX,
              segPlaneY, segPlaneZ, iChainID, 1);
    InsertArm(home, poNbr1, &newNode->myTag, -bx1, -by1, -bz1, segPlaneX,
              segPlaneY, segPlaneZ, iChainID, 1);
    InsertArm(home, newNode, &poNbr2->myTag, -bx1, -by1, -bz1, segPlaneX,
              segPlaneY, segPlaneZ, iChainID, 1);
    InsertArm(home, poNbr2, &newNode->myTag, bx1, by1, bz1, segPlaneX,
              segPlaneY, segPlaneZ, iChainID, 1);

    ChangeArmBurg(home, node, &poNbr1->myTag, 0, 0, 0, 0, 0, 0, 1, DEL_SEG_NONE,
                  true);
    ChangeArmBurg(home, poNbr1, &node->myTag, 0, 0, 0, 0, 0, 0, 1, DEL_SEG_NONE,
                  true);
    ChangeArmBurg(home, node, &poNbr2->myTag, 0, 0, 0, 0, 0, 0, 1, DEL_SEG_NONE,
                  true);
    ChangeArmBurg(home, poNbr2, &node->myTag, 0, 0, 0, 0, 0, 0, 1, DEL_SEG_NONE,
                  true);
  }
}

/*---------------------------------------------------------------------------
 *
 *      Function:       ExemptCollisionAfterDissociation
 *      Description:    For all native nodes, exempt collisions for those who
 *just dissociate onto the twin plane junjie
 *
 *-------------------------------------------------------------------------*/
void ExemptCollisionAfterDissociation(Home_t *home) {
  if (home->param->EnableTwinPlaneCrossSlip == 0) {
    return;
  }
  int i;
  Node_t *node;
  Node_t *nbrNode1, *nbrNode2, *nbrNode3;
  bool bHasOutOfPlaneSeg, bHasInPlaneSeg;
  real8 segPlaneX, segPlaneY, segPlaneZ;
  real8 A, B, C, len;
  real8 dTol = 1E-6;
  A = home->param->A;
  B = home->param->B;
  C = home->param->C;
  len = sqrt(pow(A, 2.0) + pow(B, 2.0) + pow(C, 2.0));
  segPlaneX = A / len;
  segPlaneY = B / len;
  segPlaneZ = C / len;
  for (i = 0; i < home->newNodeKeyPtr; i++) {
    node = home->nodeKeys[i];
    if (node == NULL) {
      continue;
    }
    if (node->csState == ON_CROSS_SLIP_PLANE) {
      if (node->vX == 0 && node->vY == 0 && node->vZ == 0) {
        if (node->numNbrs != 2)
          continue;
        node->flags |= NO_COLLISIONS;
        continue;
      }
    }
  }
}

/*---------------------------------------------------------------------------
 *
 *      Function:       HingeJointCanHappen
 *      Description:    Decide if a hinge-joint collisions can happen or not
 *      junjie
 *
 *-------------------------------------------------------------------------*/
bool HingeJointCanHappen(Home_t *home, Node_t *poNode1, Node_t *poNode2,
                         Node_t *poNode3) {
  if (home->param->EnableTwinPlaneCrossSlip == 0)
    return true;
  if (poNode1->csState != AT_INTERSECTION &&
      poNode1->csState != ON_CROSS_SLIP_PLANE)
    return true;
  Vector b, b2, b3;
  real8 bx2, by2, bz2;
  real8 bx3, by3, bz3;
  real8 nx, ny, nz;
  real8 dTol = 1E-6;
  GetBurgersVectorNormal(home, poNode1, poNode2, &bx2, &by2, &bz2, &nx, &ny,
                         &nz);
  GetBurgersVectorNormal(home, poNode1, poNode3, &bx3, &by3, &bz3, &nx, &ny,
                         &nz);
  b2 = Vector(bx2, by2, bz2);
  b3 = Vector(bx3, by3, bz3);
  if (ParadisCrossSlipServer::Is112Vector(b2) &&
      ParadisCrossSlipServer::Is112Vector(b3)) {
    if ((b2 + b3).Length() > dTol &&
        ParadisCrossSlipServer::Is110Vector(b2 + b3))
      return false;
  }
  return true;
}

/*---------------------------------------------------------------------------
 *
 *      Function:       SegSegCanHappen
 *      Description:    Decide if a segment-segment collisions can happen or not
 *      junjie
 *
 *-------------------------------------------------------------------------*/
bool SegSegCanHappen(Home_t *home, Node_t *poNode1, Node_t *poNode2,
                     Node_t *poNode3, Node_t *poNode4) {
  if (home->param->EnableTwinPlaneCrossSlip == 0)
    return true;
  if (poNode1->csState != ON_CROSS_SLIP_PLANE &&
      poNode1->csState != AT_INTERSECTION)
    return true; // all the 4 nodes are on the twin plane
  // make sure 1->2 3->4 are parallel to each other
  real8 dx1, dy1, dz1, dx2, dy2, dz2;
  real8 dTol = 1E-6;
  dx1 = poNode2->x - poNode1->x;
  dy1 = poNode2->y - poNode1->y;
  dz1 = poNode2->z - poNode1->z;
  dx2 = poNode4->x - poNode3->x;
  dy2 = poNode4->y - poNode3->y;
  dz2 = poNode4->z - poNode3->z;
  Vector d1, d2;
  d1 = Vector(dx1, dy1, dz1);
  d2 = Vector(dx2, dy2, dz2);
  if (d1 * d2 < 0) {
    Node_t *tmp;
    tmp = poNode3;
    poNode3 = poNode4;
    poNode4 = tmp;
  }
  // make sure the burgers vector of 1->2 and 3->4 are of [112] type
  Vector b1, b2;
  real8 bx1, by1, bz1;
  real8 bx2, by2, bz2;
  real8 nx, ny, nz;
  GetBurgersVectorNormal(home, poNode1, poNode2, &bx1, &by1, &bz1, &nx, &ny,
                         &nz);
  GetBurgersVectorNormal(home, poNode3, poNode4, &bx2, &by2, &bz2, &nx, &ny,
                         &nz);
  b1 = Vector(bx1, by1, bz1);
  b2 = Vector(bx2, by2, bz2);
  if (ParadisCrossSlipServer::Is112Vector(b1) &&
      ParadisCrossSlipServer::Is112Vector(b2)) {
    // if ((b1+b2).Length() > dTol) return false;
    if ((b1 + b2).Length() > dTol &&
        ParadisCrossSlipServer::Is110Vector(b1 + b2))
      return false;
  }
  return true;
}

/*---------------------------------------------------------------------------
 *
 *      Function:       SurfaceNodeDissociation
 *      Description:    Dissociate the intersection nodes on the surface
 *
 *      junjie
 *
 *-------------------------------------------------------------------------*/
void SurfaceNodeDissociation(Home_t *home, Node_t *node) {
  real8 lenDissociate = home->param->rmax + 1.0;
  Node_t *newNode, *poNbr1, *poNbr2;
  ;
  real8 bx, by, bz, nx, ny, nz;
  real8 lenDeduc;
  real8 pos[3];
  real8 segPlaneX, segPlaneY, segPlaneZ;
  real8 A, B, C, len;
  real8 daTemp[3];
  real8 dx, dy, dz;
  int iChainID = 0;
  A = home->param->A;
  B = home->param->B;
  C = home->param->C;
  len = sqrt(pow(A, 2.0) + pow(B, 2.0) + pow(C, 2.0));
  segPlaneX = A / len;
  segPlaneY = B / len;
  segPlaneZ = C / len;
  Node_t *poNbr = GetNeighborNode(home, node, 0);
  if (poNbr == (Node_t *)NULL) {
    Fatal("Neighbor not found at %s line %d\n", __FILE__, __LINE__);
  }
  GetBurgersVectorNormal(home, node, poNbr, &bx, &by, &bz, &nx, &ny, &nz);
  Vector ovelo(node->vX, node->vY, node->vZ);
  Vector oNormal1(segPlaneX, segPlaneY, segPlaneZ);
  Vector oNormal2(nx, ny, nz);
  Vector oSur(node->dNx, node->dNy, node->dNz);
  Vector oline = oNormal1 ^ oNormal2;
  lenDeduc = oline * ovelo;
  ovelo = ovelo - oline * lenDeduc;
  lenDeduc = ovelo * oNormal1;
  Vector oPointToLeading = ovelo - oNormal1 * lenDeduc;
  lenDeduc = oSur * oPointToLeading;
  oPointToLeading = oPointToLeading - oSur * lenDeduc;
  oPointToLeading.Normalize();
  pos[X] = node->x + lenDissociate * oPointToLeading.GetX();
  pos[Y] = node->y + lenDissociate * oPointToLeading.GetY();
  pos[Z] = node->z + lenDissociate * oPointToLeading.GetZ();
  daTemp[0] = node->x;
  daTemp[1] = node->y;
  daTemp[2] = node->z;
  poNbr->dNx = node->dNx;
  poNbr->dNy = node->dNy;
  poNbr->dNz = node->dNz;
  ChangeArmBurg(home, node, &poNbr->myTag, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1,
                DEL_SEG_NONE, true);
  ChangeArmBurg(home, poNbr, &node->myTag, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1,
                DEL_SEG_NONE, true);
  RemoveNode(home, node, 1);
  RepositionNode(home, daTemp, &(poNbr->myTag), 1);
  poNbr->constraint = SURFACE_NODE;
  poNbr->csState = ON_CROSS_SLIP_PLANE;
  newNode = GetNewNativeNode(home);
  FreeNodeArms(newNode);
  newNode->constraint = SURFACE_NODE;
  newNode->dNx = poNbr->dNx;
  newNode->dNy = poNbr->dNy;
  newNode->dNz = poNbr->dNz;
  newNode->oldx = poNbr->oldx;
  newNode->oldy = poNbr->oldy;
  newNode->oldz = poNbr->oldz;
  newNode->csState = ON_CROSS_SLIP_PLANE;
  newNode->x = pos[X];
  newNode->y = pos[Y];
  newNode->z = pos[Z];

  poNbr1 = GetNeighborNode(home, poNbr, 0);
  poNbr2 = GetNeighborNode(home, poNbr, 1);
  if (poNbr1 == (Node_t *)NULL || poNbr2 == (Node_t *)NULL) {
    Fatal("Neighbor not found at %s line %d\n", __FILE__, __LINE__);
  }

  dx = poNbr2->x - poNbr1->x;
  dy = poNbr2->y - poNbr1->y;
  dz = poNbr2->z - poNbr1->z;
  if (1) {
    AddOp(home, SPLIT_NODE, node->myTag.domainID, node->myTag.index,
          newNode->myTag.domainID, newNode->myTag.index, -1,
          -1,            /* 3rd node tag unneeded */
          0.0, 0.0, 0.0, /* burgers vector unneeded*/
          newNode->x, newNode->y, newNode->z, 0.0, 0.0,
          0.0); // normal not needed
  }
  Vector odirect(dx, dy, dz);
  if (odirect * oPointToLeading > 0) {
    GetBurgersVectorNormal(home, poNbr, poNbr2, &bx, &by, &bz, &nx, &ny, &nz);
    InsertArm(home, newNode, &poNbr2->myTag, bx, by, bz, nx, ny, nz, iChainID,
              1);
    InsertArm(home, poNbr2, &newNode->myTag, -bx, -by, -bz, nx, ny, nz,
              iChainID, 1);
    ChangeArmBurg(home, poNbr, &poNbr2->myTag, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1,
                  DEL_SEG_NONE, true);
    ChangeArmBurg(home, poNbr2, &poNbr->myTag, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1,
                  DEL_SEG_NONE, true);

  } else {
    GetBurgersVectorNormal(home, poNbr, poNbr1, &bx, &by, &bz, &nx, &ny, &nz);
    InsertArm(home, newNode, &poNbr1->myTag, bx, by, bz, nx, ny, nz, iChainID,
              1);
    InsertArm(home, poNbr1, &newNode->myTag, -bx, -by, -bz, nx, ny, nz,
              iChainID, 1);
    ChangeArmBurg(home, poNbr, &poNbr1->myTag, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1,
                  DEL_SEG_NONE, true);
    ChangeArmBurg(home, poNbr1, &poNbr->myTag, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1,
                  DEL_SEG_NONE, true);
  }
}

/*---------------------------------------------------------------------------
 *
 *      Function:       RemeshCSplane
 *      Description:    remesh the tiny segments between two multi-arm on cs
 *plane nodes
 *
 *      junjie
 *
 *-------------------------------------------------------------------------*/
void RemeshCSplane(Home_t *home) {
  if (home->param->EnableTwinPlaneCrossSlip == 0) {
    return;
  }
  Node_t *node;
  double A, B, C, len;
  A = home->param->A;
  B = home->param->B;
  C = home->param->C;
  len = sqrt(pow(A, 2.0) + pow(B, 2.0) + pow(C, 2.0));
  real8 segPlaneX = A / len;
  real8 segPlaneY = B / len;
  real8 segPlaneZ = C / len;
  real8 dTol = 1E-6;
  Node_t *poNbr1;
  bool bAllSegOnPlane;
  for (int i = 0; i < home->newNodeKeyPtr; i++) {
    if ((node = home->nodeKeys[i]) == NULL)
      continue;
    if (node->constraint == PINNED_NODE)
      continue;
    if (node->myTag.domainID != home->myDomain)
      continue;
    if (node->csState != ON_CROSS_SLIP_PLANE)
      continue;
    if (node->numNbrs <= 2)
      continue;
    bool bAllSegOnPlane = true;
    for (int k = 0; k < node->numNbrs; ++k) {
      poNbr1 = GetNeighborNode(home, node, k);
      if (poNbr1->numNbrs <= 2)
        continue;
      bAllSegOnPlane = true;
      for (int j = 0; j < poNbr1->numNbrs; ++j) {
        if (pow(node->nx[j] - segPlaneX, 2) + pow(node->ny[j] - segPlaneY, 2) +
                pow(node->nz[j] - segPlaneZ, 2) >
            dTol)
          bAllSegOnPlane = false;
      }
      if (bAllSegOnPlane)
        break;
    }

    real8 dx, dy, dz;
    dx = node->x - poNbr1->x;
    dy = node->y - poNbr1->y;
    dz = node->z - poNbr1->z;
    real8 dNewPosition[3];
    Node_t *poTargetNode = NULL;
    dNewPosition[0] = (node->x + poNbr1->x) / 2.0;
    dNewPosition[1] = (node->y + poNbr1->y) / 2.0;
    dNewPosition[2] = (node->z + poNbr1->z) / 2.0;
    if (sqrt(pow(dx, 2) + pow(dy, 2) + pow(dz, 2)) <= home->param->minSeg / 2) {
      // Merge node and poNbr1
      MergeNode(home, OPCLASS_COLLISION, node, poNbr1, dNewPosition,
                &poTargetNode, 1);
    }
  }
}