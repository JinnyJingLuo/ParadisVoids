// Ahmed M. Hussein

#include "ParadisSurface.h"
#include "list"
#include "DSMPI.h"
#include "Comm.h"
#include "float.h"
#include "Randomizer.h"
#include "QuadPatch.h"
#include "TriPatch.h"

using namespace std;
using namespace EZ;
using namespace GeometrySystem;

ParadisSurface::ParadisSurface() { Initialize(); }
ParadisSurface::ParadisSurface(const ParadisSurface &oSurface) {
  *this = oSurface;
}
ParadisSurface::~ParadisSurface() { Reset(); }
ParadisSurface &ParadisSurface::operator=(const ParadisSurface &oSurface) {
  Reset();
  m_dXMin = oSurface.m_dXMin;
  m_dXMax = oSurface.m_dXMax;
  m_dYMin = oSurface.m_dYMin;
  m_dYMax = oSurface.m_dYMax;
  m_dZMin = oSurface.m_dZMin;
  m_dZMax = oSurface.m_dZMax;
  m_dVolume = oSurface.m_dVolume;
  m_dTolerance = oSurface.m_dTolerance;
  GenerateTriangulations();
  return *this;
}
void ParadisSurface::Reset() {
  ClearTriangulations();
  list<SurfaceSegment *>::iterator liSegments;
  for (liSegments = m_lpoSurfaceSegments.begin();
       liSegments != m_lpoSurfaceSegments.end(); liSegments++) {
    if ((*liSegments) != NULL) {
      delete (*liSegments);
    }
  }
  m_lpoSurfaceSegments.clear();
}
void ParadisSurface::CheckNodes(Home_t *poHome) {
  // printf("Checked!\n");
  unsigned int iSize = (unsigned int)poHome->newNodeKeyPtr;
  unsigned int i = 0;
  Node_t *poNode = NULL;
  Point oNodePoint;
  Point oNearestSurfacePoint;
  TriPatch *poTriangle = NULL;
  Vector oSurfaceNormal;
  double daPos[3] = {0.0, 0.0, 0.0};
  ClearOpList(poHome);
  unsigned int iConstraintType;
  Vector oConstraintVector;
  double dSurfaceAnnihilationDistance = 1.0;
  bool bIsSurfaceNode = false;
  bool bIsVoidSurfaceNode = false;
  double dX = 0.0;
  double dY = 0.0;
  double dZ = 0.0;
  double temp_r = 0.0;
  for (i = 0; i < iSize; i++) {
    poNode = poHome->nodeKeys[i];
    if (poNode == NULL) {
      continue;
    }
    if (poNode->constraint == PINNED_NODE) {
      continue;
    }
    /*if(poNode->constraint > dX || poNode->constraint >dY || poNode->constraint
    >dZ   )
    {
            break;
    }*/
    dX = poNode->x;
    dY = poNode->y;
    dZ = poNode->z;
    GetInLocalCoordinates(poHome->param, dX, dY, dZ);
    oNodePoint.Set(dX, dY, dZ);
    // a node is converted to a surface node if it is outside the box or the
    // minimum distance to the surface is less than the surface annihilation
    // distance
    bIsSurfaceNode = false;
    bIsVoidSurfaceNode = false;
    // Yejun
    if (poHome->param->GeometryType == 1) {
      if (fabs(dX) > 0.5 * poHome->param->Dimensions[0] ||
          fabs(dY) > 0.5 * poHome->param->Dimensions[1] ||
          fabs(dZ) > 0.5 * poHome->param->Dimensions[2])
        bIsSurfaceNode = true;
    } else {
      bIsSurfaceNode = !IsPointInside(oNodePoint);
    }
    bIsSurfaceNode = bIsSurfaceNode || (GetLeastSurfaceDistance(oNodePoint) <
                                        dSurfaceAnnihilationDistance);
    bIsVoidSurfaceNode = (((pow(pow(poNode->x,2.0)+pow(poNode->y,2.0)+pow(poNode->z,2.0),0.5))<poHome->param->voidR+dSurfaceAnnihilationDistance) || (poNode->constraint == VOID_NODE));
    // printf("yes!\n");
    if (bIsSurfaceNode) {
      if (poHome->param->BoundaryType == RIGID_BOUNDARY)
        HandleRigidBoundary(poHome, poNode);
      else {
        GetNodeDynamicConstraint(poNode, iConstraintType, oConstraintVector);
        // the constraint vector needs to be transformed to the local
        // coordinates as well
        dX = oConstraintVector.GetX();
        dY = oConstraintVector.GetY();
        dZ = oConstraintVector.GetZ();
        GetInLocalCoordinates(poHome->param, dX, dY, dZ);
        oConstraintVector.Set(dX, dY, dZ);
        poTriangle = NULL;
        if (iConstraintType == 0) {
          // node is free to move in the space
          poTriangle = GetNearestTriangle(oNodePoint, oNearestSurfacePoint);
        } else if (iConstraintType == 1) {
          // node moves in a plane
          Plane oPlane(oConstraintVector, oNodePoint);
          poTriangle = GetNearestTriangleOnPlane(oPlane, oNodePoint,
                                                 oNearestSurfacePoint);
        } else if (iConstraintType == 2) {
          // node moves on a line
          Line oLine(oConstraintVector, oNodePoint);
          poTriangle =
              GetNearestTriangleOnLine(oLine, oNodePoint, oNearestSurfacePoint);
        } else if (iConstraintType == 3) {
          // node cannot move at all
          oNearestSurfacePoint = oNodePoint;
        }
        if (poTriangle == NULL) {
          // in case the triangle was not obtained at all, just
          // get the nearest surface point. this is not the best
          // solution in this case.
          poTriangle = GetNearestTriangle(oNodePoint, oNearestSurfacePoint);
        }

        if (poHome->param->BoundaryType == ONLY_PZ_RIGID_BOUNDARY) {
          // see if the node has escaped from the positive z surface
          if (oNodePoint.GetZ() >
              0.5 *
                  poHome->param
                      ->Dimensions[2]) //- dSurfaceAnnihilationDistance)
          {
            daPos[0] = oNearestSurfacePoint.GetX();
            daPos[1] = oNearestSurfacePoint.GetY();
            daPos[2] = oNearestSurfacePoint.GetZ();
            GetInGlobalCoordinates(poHome->param, daPos[0], daPos[1], daPos[2]);
            RepositionNode(poHome, daPos, &(poNode->myTag), 1);
            ChangeNodeType(poHome, poNode, SURFACE_NODE, 1);
            // get the surface normal in the global coordinates
            oSurfaceNormal = poTriangle->GetNormal();
            dX = oSurfaceNormal.GetX();
            dY = oSurfaceNormal.GetY();
            dZ = oSurfaceNormal.GetZ();
            GetInGlobalCoordinates(poHome->param, dX, dY, dZ);
            poNode->dNx = dX;
            poNode->dNy = dY;
            poNode->dNz = dZ;
          } else
            HandleRigidBoundary(poHome, poNode);
        } else if (poHome->param->BoundaryType == FREE_BOUNDARY) {
          daPos[0] = oNearestSurfacePoint.GetX();
          daPos[1] = oNearestSurfacePoint.GetY();
          daPos[2] = oNearestSurfacePoint.GetZ();
          GetInGlobalCoordinates(poHome->param, daPos[0], daPos[1], daPos[2]);
          RepositionNode(poHome, daPos, &(poNode->myTag), 1);
          if (poHome->param->BoundaryType == FREE_BOUNDARY) {
            ChangeNodeType(poHome, poNode, SURFACE_NODE, 1);
          }
          // get the surface normal in the global coordinates
          oSurfaceNormal = poTriangle->GetNormal();
          dX = oSurfaceNormal.GetX();
          dY = oSurfaceNormal.GetY();
          dZ = oSurfaceNormal.GetZ();
          GetInGlobalCoordinates(poHome->param, dX, dY, dZ);
          poNode->dNx = dX;
          poNode->dNy = dY;
          poNode->dNz = dZ;
        }
      }
    }
  if (bIsVoidSurfaceNode)
  {
    // printf("Inside!");
    temp_r=(pow(pow(poNode->x,2.0)+pow(poNode->y,2.0)+pow(poNode->z,2.0),0.5));
    if (temp_r<1e-6)
    {
      daPos[0] = poNode->x/(temp_r+1)*poHome->param->voidR;
      daPos[1] = poNode->y/(temp_r+1)*poHome->param->voidR;
      daPos[2] = poNode->z/(temp_r+1)*poHome->param->voidR;
    }
    else
    {
      daPos[0] = poNode->x/temp_r*poHome->param->voidR;
      daPos[1] = poNode->y/temp_r*poHome->param->voidR;
      daPos[2] = poNode->z/temp_r*poHome->param->voidR;
    }
    //GetInGlobalCoordinates(poHome->param, daPos[0], daPos[1], daPos[2]);
    RepositionNode(poHome, daPos, &(poNode->myTag), 1);
    ChangeNodeType(poHome, poNode, VOID_NODE, 1);
    dX = poNode->x;
    dY = poNode->y;
    dZ = poNode->z;
    GetInGlobalCoordinates(poHome->param, dX, dY, dZ);
    poNode->dNx = dX;
    poNode->dNy = dY;
    poNode->dNz = dZ;    
  }
  }
  CommSendRemesh(poHome);
  FixRemesh(poHome);
  DSMPI::Barrier();
}
void ParadisSurface::StoreSurfaceArms(Home_t *poHome) {
  DSMPI::Barrier();
  ClearOpList(poHome);
  unsigned int iArmsCount = 0;
  unsigned int i = 0;
  unsigned int j = 0;
  unsigned int iSize = (unsigned int)poHome->newNodeKeyPtr;
  Node_t *poNode = NULL;
  Node_t *poNeighbour = NULL;
  bool bRemovableNode = false;
  DislocationSegment oSegment;
  DislocationNode oNode1;
  DislocationNode oNode2;
  SurfaceSegment *poSurfaceSegment = NULL;
  unsigned int iRemovedArmsCount = 0;
  double dTolerance = 1.0E-3;
  for (i = 0; i < iSize; i++) {
    poNode = poHome->nodeKeys[i];
    if (poNode == NULL) {
      continue;
    }
    if (poNode->constraint >=SURFACE_NODE) {
      iArmsCount = poNode->numNbrs;
      iRemovedArmsCount = 0;
      for (j = 0; j < iArmsCount; j++) {
        poNeighbour = GetNeighborNode(poHome, poNode, j - iRemovedArmsCount);
        if (poNeighbour == NULL) {
          continue;
        }
        // don't double count surface segments
        if (fabs(poNode->burgX[j]) < dTolerance) {
          if (fabs(poNode->burgY[j]) < dTolerance) {
            if (poNode->burgZ[j] <= dTolerance) {
              continue;
            }
          } else if (poNode->burgY[j] < 0.0) {
            continue;
          }
        } else if (poNode->burgX[j] < 0.0) {
          continue;
        }
        if (poNeighbour->constraint >=SURFACE_NODE) {
          // store the segment in the surface segments list
          oSegment.SetSlipPlaneNormal(
              Vector(poNode->nx[j], poNode->ny[j], poNode->nz[j]));
          oSegment.SetBurgersVector(
              Vector(poNode->burgX[j], poNode->burgY[j], poNode->burgZ[j]));
          oNode1.Set(poNode->x, poNode->y, poNode->z);
          oNode2.Set(poNeighbour->x, poNeighbour->y, poNeighbour->z);
          oNode1.SetSurfaceNormal(
              Vector(poNode->dNx, poNode->dNy, poNode->dNz));
          oNode2.SetSurfaceNormal(
              Vector(poNeighbour->dNx, poNeighbour->dNy, poNeighbour->dNz));
          poSurfaceSegment = new SurfaceSegment;
          poSurfaceSegment->Set(oNode1, oNode2, oSegment);
          m_lpoSurfaceSegments.push_back(poSurfaceSegment);
          // and remove it from the system
          ChangeArmBurg(poHome, poNode, &poNeighbour->myTag, 0.0, 0.0, 0.0, 0.0,
                        0.0, 0.0, 1, DEL_SEG_NONE, true);
          ChangeArmBurg(poHome, poNeighbour, &poNode->myTag, 0.0, 0.0, 0.0, 0.0,
                        0.0, 0.0, 1, DEL_SEG_NONE, true);
          iRemovedArmsCount = iRemovedArmsCount + 1;
        }
      }
    }
  }
  // remove all the non surface nodes that are connected only to surface nodes
  for (i = 0; i < iSize; i++) {
    poNode = poHome->nodeKeys[i];
    if (poNode == NULL) {
      continue;
    }
    if (poNode->constraint != SURFACE_NODE) {
      // if the node is not a surface node, but connected to surface nodes only,
      // remove all of its arms and remove the node itself
      iArmsCount = poNode->numNbrs;
      bRemovableNode = true;
      for (j = 0; j < iArmsCount; j++) {
        poNeighbour = GetNeighborNode(poHome, poNode, j);
        if (poNeighbour == NULL) {
          bRemovableNode = false;
          break;
        }
        if (poNeighbour->constraint != SURFACE_NODE) {
          bRemovableNode = false;
          break;
        }
      }
      if (bRemovableNode) {
        for (j = 0; j < iArmsCount; j++) {
          // it is always zero because the all of the arms are being removed
          poNeighbour = GetNeighborNode(poHome, poNode, 0);
          if (poNeighbour == NULL) {
            continue;
          }
          if (poNeighbour->constraint >=SURFACE_NODE) {
            ChangeArmBurg(poHome, poNode, &poNeighbour->myTag, 0.0, 0.0, 0.0,
                          0.0, 0.0, 0.0, 1, DEL_SEG_NONE, true);
            ChangeArmBurg(poHome, poNeighbour, &poNode->myTag, 0.0, 0.0, 0.0,
                          0.0, 0.0, 0.0, 1, DEL_SEG_NONE, true);
          }
        }
      }
    }
  }

  // get rid of any nodes that have no arms
  CommSendRemesh(poHome);
  FixRemesh(poHome);
  RemoveOrphanedNodes(poHome);
  CommSendRemesh(poHome);
  FixRemesh(poHome);
  DSMPI::Barrier();
}
void ParadisSurface::StoreSurfaceNodesMotion(Home_t *poHome) {
  unsigned int i = 0;
  unsigned int iSize = (unsigned int)poHome->newNodeKeyPtr;
  Node_t *poNode = NULL;
  double dX = 0.0;
  double dY = 0.0;
  double dZ = 0.0;
  double dTolerance = 1.0E-6;
  DislocationSegment oSegment;
  DislocationNode oNode1;
  DislocationNode oNode2;
  SurfaceSegment *poSurfaceSegment = NULL;
  for (i = 0; i < iSize; i++) {
    poNode = poHome->nodeKeys[i];
    if (poNode == NULL) {
      continue;
    }
    // the node has to be a surface node with 1 arm
    if ((poNode->constraint >=SURFACE_NODE) && (poNode->numNbrs == 1)) {
      dX = poNode->oldx - poNode->x;
      dY = poNode->oldy - poNode->y;
      dZ = poNode->oldz - poNode->z;
      // and the motion has to be for a nonzero distance
      if ((dX * dX + dY * dY + dZ * dZ) > dTolerance) {
        oSegment.SetSlipPlaneNormal(
            Vector(poNode->nx[0], poNode->ny[0], poNode->nz[0]));
        oSegment.SetBurgersVector(
            Vector(poNode->burgX[0], poNode->burgY[0], poNode->burgZ[0]));
        oNode1.Set(poNode->oldx, poNode->oldy, poNode->oldz);
        oNode2.Set(poNode->x, poNode->y, poNode->z);
        oNode1.SetSurfaceNormal(Vector(poNode->dNx, poNode->dNy, poNode->dNz));
        oNode2.SetSurfaceNormal(Vector(poNode->dNx, poNode->dNy, poNode->dNz));
        poSurfaceSegment = new SurfaceSegment;
        poSurfaceSegment->Set(oNode1, oNode2, oSegment);
        m_lpoSurfaceSegments.push_back(poSurfaceSegment);
      }
    }
  }
}
void ParadisSurface::WriteSurfaceSegments(Home_t *poHome) {
  // reduce the segments in the surface segments list
  ReduceSurfaceSemgnts();
  // open the surface segments file
  char cWrite[64];
  sprintf(cWrite, "surface_segments_%d.txt", poHome->myDomain);
  FILE *fpFile = fopen(cWrite, "a");
  // write any segments in the list and delete them as we write
  list<SurfaceSegment *>::iterator liSegments;
  const Vector *poBurgersVector = NULL;
  const Vector *poSlipPlaneNormal = NULL;
  const DislocationNode *poFirstNode = NULL;
  const DislocationNode *poSecondNode = NULL;
  Vector oFirstSurfaceNormal;
  Vector oSecondSurfaceNormal;
  for (liSegments = m_lpoSurfaceSegments.begin();
       liSegments != m_lpoSurfaceSegments.end(); liSegments++) {
    if ((*liSegments) != NULL) {
      // the file writing format
      // time step, processor id, node1 tag, node2 tag, node1 position, node2
      // position node1 surface normal, node2 surface normal, segment burgers
      // vector, segment slip plane
      poBurgersVector = (*liSegments)->GetBurgersVector();
      poSlipPlaneNormal = (*liSegments)->GetSlipPlaneNormal();
      poFirstNode = (*liSegments)->GetFirstNode();
      poSecondNode = (*liSegments)->GetSecondNode();
      oFirstSurfaceNormal = poFirstNode->GetSurfaceNormal();
      oSecondSurfaceNormal = poSecondNode->GetSurfaceNormal();
      fprintf(fpFile,
              "%d,%d : (0,0) -> (0,0) : (%lf,%lf,%lf) / (%lf,%lf,%lf)\n",
              poHome->cycle, poHome->myDomain, poSlipPlaneNormal->GetX(),
              poSlipPlaneNormal->GetY(), poSlipPlaneNormal->GetZ(),
              poBurgersVector->GetX(), poBurgersVector->GetY(),
              poBurgersVector->GetZ());
      fprintf(
          fpFile,
          "(%lf,%lf,%lf) -> (%lf,%lf,%lf) : (%lf,%lf,%lf) -> (%lf,%lf,%lf)\n",
          poFirstNode->GetX(), poFirstNode->GetY(), poFirstNode->GetZ(),
          poSecondNode->GetX(), poSecondNode->GetY(), poSecondNode->GetZ(),
          oFirstSurfaceNormal.GetX(), oFirstSurfaceNormal.GetY(),
          oFirstSurfaceNormal.GetZ(), oSecondSurfaceNormal.GetX(),
          oSecondSurfaceNormal.GetY(), oSecondSurfaceNormal.GetZ());
      delete (*liSegments);
    }
  }
  m_lpoSurfaceSegments.clear();
  fclose(fpFile);
}
void ParadisSurface::ReduceSurfaceSemgnts() {
  unsigned int iInitialSegmentsCount =
      (unsigned int)m_lpoSurfaceSegments.size();
  list<SurfaceSegment *>::iterator liOuterSegments;
  list<SurfaceSegment *>::iterator liInnerSegments;
  double dTolerance = 1.0E-6;
  liOuterSegments = m_lpoSurfaceSegments.begin();
  while (liOuterSegments != m_lpoSurfaceSegments.end()) {
    if ((*liOuterSegments)->GetLength() < dTolerance) {
      delete (*liOuterSegments);
      liOuterSegments = m_lpoSurfaceSegments.erase(liOuterSegments);
    } else {
      liInnerSegments = liOuterSegments;
      liInnerSegments++;
      while (liInnerSegments != m_lpoSurfaceSegments.end()) {
        if ((*liOuterSegments)->IsPostMergeable((*liInnerSegments))) {
          (*liOuterSegments)->PostMerge((*liInnerSegments));
          delete (*liInnerSegments);
          liInnerSegments = m_lpoSurfaceSegments.erase(liInnerSegments);
        } else {
          liInnerSegments++;
        }
      }
      liOuterSegments++;
    }
  }
  unsigned int iFinalSegmentsCount = (unsigned int)m_lpoSurfaceSegments.size();
}
void ParadisSurface::Initialize() {
  m_lpoTriangles.clear();
  m_lpoPoints.clear();
  m_lpoSurfaceSegments.clear();
  m_dXMin = 0.0;
  m_dXMax = 0.0;
  m_dYMin = 0.0;
  m_dYMax = 0.0;
  m_dZMin = 0.0;
  m_dZMax = 0.0;
  m_dVolume = 0.0;
  m_dTolerance = 0.0;
}
void ParadisSurface::GetNodeDynamicConstraint(Node_t *poNode,
                                              unsigned int &iConstraintType,
                                              Vector &oConstraintVector) {
  oConstraintVector.Set(0.0, 0.0, 0.0);
  iConstraintType = 0;
  double dTolerance = 1.0E-6;
  unsigned int i = 0;
  // get the node dynamic constraint
  if (poNode->numNbrs == 0) {
    oConstraintVector.Set(0.0, 0.0, 0.0);
    iConstraintType = 0;
  } else if (poNode->numNbrs == 1) {
    oConstraintVector.Set(poNode->nx[0], poNode->ny[0], poNode->nz[0]);
    iConstraintType = 1;
  } else if (poNode->numNbrs == 2) {
    Vector oNormal1(poNode->nx[0], poNode->ny[0], poNode->nz[0]);
    Vector oNormal2(poNode->nx[1], poNode->ny[1], poNode->nz[1]);
    oNormal1.Normalize();
    oNormal2.Normalize();
    Vector oCrossProduct = oNormal1 ^ oNormal2;
    if (oCrossProduct.Length() < dTolerance) {
      iConstraintType = 1;
      oConstraintVector = oNormal1;
    } else {
      iConstraintType = 2;
      oConstraintVector = oCrossProduct;
    }
  } else {
    list<Vector> loNormals;
    list<Vector>::iterator liNormals;
    loNormals.push_back(Vector(poNode->nx[0], poNode->ny[0], poNode->nz[0]));
    loNormals.front().Normalize();
    Vector oTempNormal;
    real8 oDotProduct;
    bool bAddNormal = false;
    for (i = 1; i < poNode->numNbrs; i++) {
      oTempNormal.Set(poNode->nx[i], poNode->ny[i], poNode->nz[i]);
      oTempNormal.Normalize();
      bAddNormal = true;
      for (liNormals = loNormals.begin(); liNormals != loNormals.end();
           liNormals++) // junjie: found a problem here!
      {
        oDotProduct = oTempNormal * (*liNormals);
        if (fabs(1.0 - fabs(oDotProduct)) < dTolerance) {
          bAddNormal = false;
          break;
        }
      }
      if (bAddNormal) {
        loNormals.push_back(oTempNormal);
      }
    }
    unsigned int iDynamicConstraintsCount = (unsigned int)loNormals.size();
    if (iDynamicConstraintsCount == 1) {
      iConstraintType = 1;
      oConstraintVector.Set(poNode->nx[0], poNode->ny[0], poNode->nz[0]);
    } else if (iDynamicConstraintsCount == 2) {
      iConstraintType = 2;
      oConstraintVector = loNormals.front() ^ loNormals.back();
    } else {
      iConstraintType = 3;
      oConstraintVector.Set(0.0, 0.0, 0.0);
    }
  }
  oConstraintVector.Normalize();
}
TriPatch *ParadisSurface::GetNearestTriangle(const Point &oPoint,
                                             Point &oNearestPoint) const {
  list<TriPatch *>::const_iterator liTriangles;
  double dTemp = 0.0;
  double dMinDistance = DBL_MAX;
  TriPatch *poTriangle = NULL;
  Point oTempPoint;
  for (liTriangles = m_lpoTriangles.begin();
       liTriangles != m_lpoTriangles.end(); liTriangles++) {
    oTempPoint = (*liTriangles)->GetNearestPoint(oPoint, dTemp);
    if (dTemp < dMinDistance) {
      dMinDistance = dTemp;
      poTriangle = *liTriangles;
      oNearestPoint = oTempPoint;
    }
  }
  return poTriangle;
}
TriPatch *ParadisSurface::GetNearestTriangleOnPlane(
    const Plane &oPlane, const Point &oPoint, Point &oNearestPoint) const {
  list<TriPatch *>::const_iterator liTriangles;
  double dTemp = 0.0;
  double dMinDistance = DBL_MAX;
  TriPatch *poTriangle = NULL;
  Point oTempPoint;
  for (liTriangles = m_lpoTriangles.begin();
       liTriangles != m_lpoTriangles.end(); liTriangles++) {
    if ((*liTriangles)
            ->GetNearestPointOnPlane(oPoint, oPlane, oTempPoint, dTemp)) {
      if (dTemp < dMinDistance) {
        dMinDistance = dTemp;
        poTriangle = *liTriangles;
        oNearestPoint = oTempPoint;
      }
    }
  }
  return poTriangle;
}
TriPatch *ParadisSurface::GetNearestTriangleOnLine(const Line &oLine,
                                                   const Point &oPoint,
                                                   Point &oNearestPoint) const {
  list<TriPatch *>::const_iterator liTriangles;
  double dTemp = 0.0;
  double dMinDistance = DBL_MAX;
  TriPatch *poTriangle = NULL;
  Point oTempPoint;
  for (liTriangles = m_lpoTriangles.begin();
       liTriangles != m_lpoTriangles.end(); liTriangles++) {
    if ((*liTriangles)
            ->GetNearestPointOnLine(oPoint, oLine, oTempPoint, dTemp)) {
      if (dTemp < dMinDistance) {
        dMinDistance = dTemp;
        poTriangle = *liTriangles;
        oNearestPoint = oTempPoint;
      }
    }
  }
  return poTriangle;
}
bool ParadisSurface::IsPointInside(const Point &oPoint) const {
  unsigned int iTestPointsCount = 5;
  double dAmplificationFactor = 1.5;
  unsigned int i = 0;
  unsigned int j = 0;
  Point oTestPoint;
  double dMaxDimension =
      max((m_dXMax - m_dXMin), max((m_dYMax - m_dYMin), (m_dZMax - m_dZMin)));
  double dXMean = 0.5 * (m_dXMax + m_dXMin);
  double dYMean = 0.5 * (m_dYMax + m_dYMin);
  double dZMean = 0.5 * (m_dZMax + m_dZMin);
  Point oCenter(dXMean, dYMean, dZMean);
  double dRadius = dAmplificationFactor * dMaxDimension;
  double dTheta = 0.0;
  double dPhi = 0.0;
  list<TriPatch *>::const_iterator liTriangles;
  double dIntersectionParameter = 0.0;
  Point oIntersectionPoint;
  unsigned int iCount = 0;
  unsigned int iIsPointInsideCount = 0;
  for (i = 0; i < iTestPointsCount; i++) {
    dPhi = Randomizer::Random(0, PI);
    dTheta = Randomizer::Random(0, 2 * PI);
    oTestPoint.Set(dRadius * sin(dPhi) * cos(dTheta),
                   dRadius * sin(dPhi) * sin(dTheta), dRadius * cos(dPhi));
    oTestPoint = oTestPoint + oCenter;
    iCount = 0;
    for (liTriangles = m_lpoTriangles.begin();
         liTriangles != m_lpoTriangles.end(); liTriangles++) {
      if ((*liTriangles)
              ->GetLineIntersection(oPoint, oTestPoint, oIntersectionPoint,
                                    dIntersectionParameter)) {
        iCount = iCount + 1;
      }
    }
    if (iCount % 2 != 0) {
      iIsPointInsideCount = iIsPointInsideCount + 1;
    }
  }
  if (2 * iIsPointInsideCount > iTestPointsCount) {
    return true;
  }
  return false;
}
bool ParadisSurface::IsPointOnSurface(const Point &oPoint) const {
  list<TriPatch *>::const_iterator liTriangles;
  for (liTriangles = m_lpoTriangles.begin();
       liTriangles != m_lpoTriangles.end(); liTriangles++) {
    if ((*liTriangles)->IsPointInTriangle(oPoint)) {
      return true;
    }
  }
  return false;
}
double ParadisSurface::GetLeastSurfaceDistance(const Point &oPoint) const {
  list<TriPatch *>::const_iterator liTriangles;
  double dMin = DBL_MAX;
  double dTemp = 0.0;
  for (liTriangles = m_lpoTriangles.begin();
       liTriangles != m_lpoTriangles.end(); liTriangles++) {
    dTemp = (*liTriangles)->GetNormalDistance(oPoint);
    if (dTemp < dMin) {
      dMin = dTemp;
    }
  }
  return dMin;
}
void ParadisSurface::ClearTriangulations() {
  list<TriPatch *>::iterator liTriangles;
  for (liTriangles = m_lpoTriangles.begin();
       liTriangles != m_lpoTriangles.end(); liTriangles++) {
    if ((*liTriangles) != NULL) {
      delete (*liTriangles);
    }
  }
  m_lpoTriangles.clear();

  list<GenericNode *>::iterator liPoints;
  for (liPoints = m_lpoPoints.begin(); liPoints != m_lpoPoints.end();
       liPoints++) {
    if ((*liPoints) != NULL) {
      delete (*liPoints);
    }
  }
  m_lpoPoints.clear();
}
void ParadisSurface::GenerateSurfaceTriangulationFromFEMNodes(
    const vector<FEMNode *> &vpoNodes,
    const vector<FEMElement *> &vpoElements) {
  ClearTriangulations();
  unsigned int i = 0;
  unsigned int iNodesCount = (unsigned int)vpoNodes.size();
  unsigned int iElementsCount = (unsigned int)vpoElements.size();
  vector<GenericNode *> vpoPoints;
  vector<bool> vbIsUsed;
  vpoPoints.reserve(iNodesCount +
                    2 * iElementsCount); // the 2 factor is just an estimate
  vbIsUsed.reserve(iNodesCount + 2 * iElementsCount);
  GenericNode *poNode = NULL;
  for (i = 0; i < iNodesCount; i++) {
    poNode = new GenericNode;
    poNode->Set(vpoNodes[i]->GetX(), vpoNodes[i]->GetY(), vpoNodes[i]->GetZ());
    poNode->SetID(vpoNodes[i]->GetID());
    vpoPoints.push_back(poNode);
    vbIsUsed.push_back(false);
  }

  vector<vector<unsigned int>> vviNodesIndices;
  vector<GenericNode *> vpoMidPoints;
  unsigned int j = 0;
  unsigned int k = 0;
  unsigned int iPatchesCount = 0;
  vector<GenericNode *> vpoTempPoints;
  vector<QuadPatch> voPatches;
  voPatches.reserve(2 * iElementsCount);

  // loop over the actual elements and get the quad patches
  for (i = 0; i < iElementsCount; i++) {
    vpoElements[i]->GetGeometry()->GenerateSurfacePatches(vviNodesIndices,
                                                          vpoMidPoints);
    iPatchesCount = vviNodesIndices.size();
    for (j = 0; j < iPatchesCount; j++) {
      if (vviNodesIndices[j].size() != PointsPerQuadPatch - 1 ||
          vpoMidPoints[j] == NULL) {
        continue;
      }
      vpoPoints.push_back(vpoMidPoints[j]);
      vbIsUsed.push_back(false);
      vpoTempPoints.resize(PointsPerQuadPatch);
      for (k = 0; k < PointsPerQuadPatch - 1; k++) {
        vpoTempPoints[k] = vpoPoints[vviNodesIndices[j].at(k) - 1];
        vbIsUsed[vviNodesIndices[j].at(k) - 1] = true;
      }
      vpoTempPoints[PointsPerQuadPatch - 1] = vpoPoints.back();
      vbIsUsed.back() = true;
      voPatches.push_back(QuadPatch(vpoTempPoints));
    }
  }

  // now all the used points and the quad patches are ready
  unsigned int iSize = vpoPoints.size();
  for (i = 0; i < iSize; i++) {
    if (vbIsUsed[i]) {
      m_lpoPoints.push_back(vpoPoints[i]);
    } else {
      delete vpoPoints[i];
    }
  }
  vpoPoints.clear();

  list<TriPatch *>::iterator liSurfaceTriangles;

  iSize = voPatches.size();
  vector<Patch *> vpoTempTriangles;
  unsigned int iTrianglesCount = 0;
  for (i = 0; i < iSize; i++) {
    vpoTempTriangles = voPatches[i].GenerateTriangulation();
    iTrianglesCount = vpoTempTriangles.size();
    for (j = 0; j < iTrianglesCount; j++) {
      m_lpoTriangles.push_back((TriPatch *)vpoTempTriangles[j]);
    }
  }
  unsigned int iIndex = 0;
  list<GenericNode *>::iterator liPoints;
  for (liPoints = m_lpoPoints.begin(); liPoints != m_lpoPoints.end();
       liPoints++) {
    (*liPoints)->SetID(iIndex);
    iIndex = iIndex + 1;
  }
}
void ParadisSurface::HandleRigidBoundary(Home_t *poHome, Node_t *poNode) {
  double daPos[3] = {0.0, 0.0, 0.0};
  daPos[0] = poNode->oldx;
  daPos[1] = poNode->oldy;
  daPos[2] = poNode->oldz;
  RepositionNode(poHome, daPos, &(poNode->myTag), 1);
}
