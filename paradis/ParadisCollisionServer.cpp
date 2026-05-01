#include "ParadisCollisionServer.h"
#include "Home.h"
#include "TwinPlaneCrossSlip.h"

AnnihilationSegment::AnnihilationSegment() { Initialize(); }
AnnihilationSegment::AnnihilationSegment(const AnnihilationSegment &oSegment) {
  *this = oSegment;
}
AnnihilationSegment::~AnnihilationSegment() { Reset(); }
AnnihilationSegment &AnnihilationSegment::
operator=(const AnnihilationSegment &oSegment) {
  m_oStartPoint = oSegment.m_oStartPoint;
  m_oEndPoint = oSegment.m_oEndPoint;
  m_oBurgersVector = oSegment.m_oBurgersVector;
  m_iTimeStep = oSegment.m_iTimeStep;
  return *this;
}
void AnnihilationSegment::Reset() { Initialize(); }
void AnnihilationSegment::SetPoints(const Point &oStart, const Point &oEnd) {
  m_oStartPoint = oStart;
  m_oEndPoint = oEnd;
}
void AnnihilationSegment::SetBurgersVector(const Vector &oBurgersVector) {
  m_oBurgersVector = oBurgersVector;
}
void AnnihilationSegment::SetTimeStep(const unsigned int &iTimeStep) {
  m_iTimeStep = iTimeStep;
}
void AnnihilationSegment::Write(FILE *fpFile) const {
  fprintf(fpFile, "%d : %s -> %s @ %s\n", m_iTimeStep,
          m_oStartPoint.ToString().c_str(), m_oEndPoint.ToString().c_str(),
          m_oBurgersVector.ToString().c_str());
}
void AnnihilationSegment::Initialize() {
  m_oStartPoint.Reset();
  m_oEndPoint.Reset();
  m_oBurgersVector.Reset();
  m_iTimeStep = 0;
}

ParadisCollisionServer
    *ParadisCollisionServer::m_poParadisCollisionServerInstance = NULL;
ParadisCollisionServer *ParadisCollisionServer::CreateInstance() {
  if (m_poParadisCollisionServerInstance == NULL) {
    m_poParadisCollisionServerInstance = new ParadisCollisionServer;
  }
  return m_poParadisCollisionServerInstance;
}
ParadisCollisionServer::ParadisCollisionServer() { Initialize(); }
ParadisCollisionServer::~ParadisCollisionServer() { Reset(); }
void ParadisCollisionServer::Reset() {
  list<AnnihilationSegment *>::iterator liSegments;
  for (liSegments = m_lpoAnnihilatedSegments.begin();
       liSegments != m_lpoAnnihilatedSegments.end(); liSegments++) {
    if ((*liSegments) != NULL)
      delete (*liSegments);
  }
  m_lpoAnnihilatedSegments.clear();
  m_sAnnihilatedSegmentsFileName = "";
  if (m_fpAnnihilatedSegments != NULL) {
    fclose(m_fpAnnihilatedSegments);
  }
}
void ParadisCollisionServer::Set(Home_t *poHome) {
  Reset();
  char cWrite[64];
  sprintf(cWrite, "annihilated_segments_%d.txt", poHome->myDomain);
  m_sAnnihilatedSegmentsFileName = string(cWrite);
  m_fpAnnihilatedSegments = fopen(cWrite, "a");
}
void ParadisCollisionServer::HandleCollisions(Home_t *poHome) {
  if (poHome->param->AnnihilateSegments != 0) {
    HandleAnnihilations(poHome);
    if (poHome->cycle % poHome->param->savefreq == 0) {
      WriteAnnihilatedSegments();
    }
  }
  HandleSegmentSegmentCollisions(poHome);
  HandleHingeJointCollisions(poHome);
  HandleTriangularLoops(poHome);
}
void ParadisCollisionServer::Initialize() {
  m_lpoAnnihilatedSegments.clear();
  m_sAnnihilatedSegmentsFileName = "";
  m_fpAnnihilatedSegments = NULL;
}
bool ParadisCollisionServer::IsNodeCollisionValid(Node_t *poNode) const {
  if (poNode == NULL) {
    return false;
  }
  if (poNode->flags & NO_COLLISIONS) {
    return false;
  }
  return true;
}
bool ParadisCollisionServer::AreSegmentsCloseAndApproaching(
    Home_t *poHome, Node_t *poNode1, Node_t *poNode2, Node_t *poNode3,
    Node_t *poNode4, Point &oNearPoint1, Point &oNearPoint2) const {
  oNearPoint1.Set(0.0, 0.0, 0.0);
  oNearPoint2.Set(0.0, 0.0, 0.0);
  double dLocalCoordinate1 = 0.0;
  double dLocalCoordinate2 = 0.0;
  double dToleranceSquared = 1.0E-6;
  Point oNode1(poNode1->x, poNode1->y, poNode1->z);
  Point oNode2 = GetNearestImage(poHome->param, poNode1, poNode2);
  Point oNode3 = GetNearestImage(poHome->param, poNode1, poNode3);
  Point oNode4 = GetNearestImage(poHome->param, poNode1, poNode4);

  Vector oM(oNode1, oNode2);
  Vector oN(oNode3, oNode4);
  // extremely short segments cannot collide
  double dA = oM.LengthSquared();
  if (dA < dToleranceSquared) {
    return false;
  }
  double dB = oN.LengthSquared();
  if (dB < dToleranceSquared) {
    return false;
  }
  double dC = oM * oN;
  Vector oV13(oNode1, oNode3);
  double dD = oM * oV13;
  double dE = -(oN * oV13);
  double dF = dA * dB - dC * dC;
  // get the 2 nearest points
  if (fabs(dF) < dToleranceSquared) {
    // the segments are parallel
    Vector oV14(oNode1, oNode4);
    Vector oV23(oNode2, oNode3);
    Vector oV24(oNode2, oNode4);
    dA = oV13.LengthSquared();
    dB = oV14.LengthSquared();
    dC = oV23.LengthSquared();
    dD = oV24.LengthSquared();
    dE = dA;
    dLocalCoordinate1 = 0.0;
    dLocalCoordinate2 = 0.0;
    if (dB < dE) {
      dE = dB;
      dLocalCoordinate1 = 0.0;
      dLocalCoordinate2 = 1.0;
    }
    if (dC < dE) {
      dE = dC;
      dLocalCoordinate1 = 1.0;
      dLocalCoordinate2 = 0.0;
    }
    if (dD < dE) {
      dE = dD;
      dLocalCoordinate1 = 1.0;
      dLocalCoordinate2 = 1.0;
    }
  } else {
    dLocalCoordinate1 = (dC * dE + dD * dB) / dF;
    dLocalCoordinate2 = (dE + dC * dLocalCoordinate1) / dB;
  }

  if (dLocalCoordinate1 < 0.0) {
    dLocalCoordinate1 = 0.0;
  }
  if (dLocalCoordinate1 > 1.0) {
    dLocalCoordinate1 = 1.0;
  }

  if (dLocalCoordinate2 < 0.0) {
    dLocalCoordinate2 = 0.0;
  }
  if (dLocalCoordinate2 > 1.0) {
    dLocalCoordinate2 = 1.0;
  }

  // now check if the segments are close
  oNearPoint1 = oNode1 + oM * dLocalCoordinate1;
  oNearPoint2 = oNode3 + oN * dLocalCoordinate2;
  Vector oTemp = oNearPoint1 - oNearPoint2;
  dA = oTemp.GetX();
  dB = oTemp.GetY();
  dC = oTemp.GetZ();
  // do NOT flip the node back to its original image
  // oNearPoint2 = GetNearestImage(poHome->param,poNode3,oNearPoint2);
  double dADot = poNode1->vX - poNode3->vX +
                 dLocalCoordinate1 * (poNode2->vX - poNode1->vX) -
                 dLocalCoordinate2 * (poNode4->vX - poNode3->vX);
  double dBDot = poNode1->vY - poNode3->vY +
                 dLocalCoordinate1 * (poNode2->vY - poNode1->vY) -
                 dLocalCoordinate2 * (poNode4->vY - poNode3->vY);
  double dCDot = poNode1->vZ - poNode3->vZ +
                 dLocalCoordinate1 * (poNode2->vZ - poNode1->vZ) -
                 dLocalCoordinate2 * (poNode4->vZ - poNode3->vZ);
  dD = dA * dA + dB * dB + dC * dC;
  double dDDot = 2.0 * (dA * dADot + dB * dBDot + dC * dCDot);
  double dAbsoluteMinimum = 1.0E-2;
  if (dD < (poHome->param->rann * poHome->param->rann)) {
    if (dDDot < 0.0) {
      return true;
    }
  }
  if (dD < dAbsoluteMinimum) {
    return true;
  }
  return false;
}
bool ParadisCollisionServer::AreNodesCloseAndApproaching(
    Home_t *poHome, Node_t *poNode1, Node_t *poNode2) const {
  Point oPoint1(poNode1->x, poNode1->y, poNode1->z);
  Point oPoint2 = GetNearestImage(poHome->param, poNode1, poNode2);
  double dDistance = oPoint1.Distance(oPoint2);
  double dAbsoluteMinimum = 1.0E-10;
  if (dDistance < dAbsoluteMinimum) {
    return true;
  }
  if (dDistance > poHome->param->rann) {
    return false;
  }
  Vector oV1(poNode1->vX, poNode1->vY, poNode1->vZ);
  Vector oV2(poNode2->vX, poNode2->vY, poNode2->vZ);
  Vector oRelativeVelocity = oV2 - oV1;
  Vector oLine(oPoint1, oPoint2);
  oLine.Normalize();
  if (oLine * oRelativeVelocity < 0.0) {
    return true;
  }
  return false;
}
bool ParadisCollisionServer::IsTriangleCollapsing(Home_t *poHome,
                                                  Node_t *poNode1,
                                                  Node_t *poNode2,
                                                  Node_t *poNode3) const {
  Point oNode1(poNode1->x, poNode1->y, poNode1->z);
  Point oNode2 = GetNearestImage(poHome->param, poNode1, poNode2);
  Point oNode3 = GetNearestImage(poHome->param, poNode1, poNode3);
  Vector oV1(oNode1, oNode2);
  Vector oV2(oNode1, oNode3);
  double dL1 = oV1.Length();
  double dL2 = oV2.Length();
  oV1.Normalize();
  oV2.Normalize();
  double dCosTheta = oV1 * oV2;
  if (poNode1->csState == AT_INTERSECTION) {
    if (dCosTheta < 0.98) // the angle is too large for merging
    {
      return false;
    }
    Vector oV1Dot(poNode2->vX - poNode1->vX, poNode2->vY - poNode1->vY,
                  poNode2->vZ - poNode1->vZ);
    Vector oV2Dot(poNode3->vX - poNode1->vX, poNode3->vY - poNode1->vY,
                  poNode3->vZ - poNode1->vZ);
    oV1Dot = oV1Dot * (1.0 / dL1);
    oV2Dot = oV2Dot * (1.0 / dL2);
    double dR = (oV1Dot * oV2 + oV1 * oV2Dot) -
                dCosTheta * (oV1Dot * oV1 + oV2Dot * oV2);
    if (dR > 0) // the angle is moderate but decreasing
    {
      return true;
    }
    return false;
  } else {
    if (dCosTheta < 0.85) // the angle is too large for merging
    {
      return false;
    }
    if (dCosTheta > 0.95) // the angle is too small, merging is allowed
    {
      return true;
    }
    Vector oV1Dot(poNode2->vX - poNode1->vX, poNode2->vY - poNode1->vY,
                  poNode2->vZ - poNode1->vZ);
    Vector oV2Dot(poNode3->vX - poNode1->vX, poNode3->vY - poNode1->vY,
                  poNode3->vZ - poNode1->vZ);
    oV1Dot = oV1Dot * (1.0 / dL1);
    oV2Dot = oV2Dot * (1.0 / dL2);
    double dR = (oV1Dot * oV2 + oV1 * oV2Dot) -
                dCosTheta * (oV1Dot * oV1 + oV2Dot * oV2);
    if (dR > 0) // the angle is moderate but decreasing
    {
      return true;
    }
    return false;
  }
}
bool ParadisCollisionServer::GetPlanePlaneCollisionPoint(
    Point *poNode1, const Vector &oNormal1, Point *poNode2,
    const Vector &oNormal2, Point &oCollisionPoint) const {
  oCollisionPoint.Set(0.0, 0.0, 0.0);
  // if the 2 points lie in the same plane, the collision point is basically the
  // mid point of the line connecting the two nodes
  double dTolerance = 1.0E-6;
  Vector oRelativePosition(*poNode1, *poNode2);
  Vector oMidPoint = (*poNode1 + *poNode2) * 0.5;
  if ((oNormal1 ^ oNormal2).Length() < dTolerance) {
    // first, make sure that it is the same plane, not just 2 parallel planes
    if (fabs(oNormal1 * oRelativePosition) > dTolerance) {
      return false;
    }
    oCollisionPoint = oMidPoint;
  } else {
    Matrix oSystem(5, 5);
    oSystem.Set(1, 1, 1.0);
    oSystem.Set(1, 2, 0.0);
    oSystem.Set(1, 3, 0.0);
    oSystem.Set(1, 4, oNormal1.GetX());
    oSystem.Set(1, 5, oNormal2.GetX());

    oSystem.Set(2, 1, 0.0);
    oSystem.Set(2, 2, 1.0);
    oSystem.Set(2, 3, 0.0);
    oSystem.Set(2, 4, oNormal1.GetY());
    oSystem.Set(2, 5, oNormal2.GetY());

    oSystem.Set(3, 1, 0.0);
    oSystem.Set(3, 2, 0.0);
    oSystem.Set(3, 3, 1.0);
    oSystem.Set(3, 4, oNormal1.GetZ());
    oSystem.Set(3, 5, oNormal2.GetZ());

    oSystem.Set(4, 1, oNormal1.GetX());
    oSystem.Set(4, 2, oNormal1.GetY());
    oSystem.Set(4, 3, oNormal1.GetZ());
    oSystem.Set(4, 4, 0.0);
    oSystem.Set(4, 5, 0.0);

    oSystem.Set(5, 1, oNormal2.GetX());
    oSystem.Set(5, 2, oNormal2.GetY());
    oSystem.Set(5, 3, oNormal2.GetZ());
    oSystem.Set(5, 4, 0.0);
    oSystem.Set(5, 5, 0.0);

    Matrix oRHS(5, 1);
    oRHS.Set(1, 1, oMidPoint.GetX());
    oRHS.Set(2, 1, oMidPoint.GetY());
    oRHS.Set(3, 1, oMidPoint.GetZ());
    oRHS.Set(4, 1, oNormal1 * Vector(poNode1->GetX(), poNode1->GetY(),
                                     poNode1->GetZ()));
    oRHS.Set(5, 1, oNormal2 * Vector(poNode2->GetX(), poNode2->GetY(),
                                     poNode2->GetZ()));

    Matrix oX = oSystem.Solve(oRHS);
    oCollisionPoint.SetX(oX.Get(1, 1));
    oCollisionPoint.SetY(oX.Get(2, 1));
    oCollisionPoint.SetZ(oX.Get(3, 1));
  }
  return true;
}
bool ParadisCollisionServer::GetLinePlaneCollisionPoint(
    Point *poNode1, const Vector &oDirection, Point *poNode2,
    const Vector &oNormal, Point &oCollisionPoint) const {
  oCollisionPoint.Set(0.0, 0.0, 0.0);
  double dTolerance = 1.0E-6;
  Vector oV(*poNode1, *poNode2);
  double dNumerator = oV * oNormal;
  if (fabs(dNumerator) < dTolerance) {
    oCollisionPoint = *poNode1;
    return true;
  }
  double dDenominator = oDirection * oNormal;
  if (fabs(dDenominator) < dTolerance) {
    // line is parallel to the plane (normal to the plane's normal)
    return false;
  }
  double dEta = dNumerator / dDenominator;
  oCollisionPoint.SetX(poNode1->GetX() + dEta * oDirection.GetX());
  oCollisionPoint.SetY(poNode1->GetY() + dEta * oDirection.GetY());
  oCollisionPoint.SetZ(poNode1->GetZ() + dEta * oDirection.GetZ());
  return true;
}
bool ParadisCollisionServer::GetLineLineCollisionPoint(
    Point *poNode1, const Vector &oDirection1, Point *poNode2,
    const Vector &oDirection2, Point &oCollisionPoint) const {
  // Can only handle the case when 2 points are on the same line(either
  // oDirection1 or oDirection2)
  oCollisionPoint.Set(0.0, 0.0, 0.0);
  oCollisionPoint.SetX((poNode1->GetX() + poNode2->GetX()) / 2);
  oCollisionPoint.SetY((poNode1->GetY() + poNode2->GetY()) / 2);
  oCollisionPoint.SetZ((poNode1->GetZ() + poNode2->GetZ()) / 2);
  double dTolerance = 1.0E-6;
  Vector oV(*poNode1, *poNode2);
  oV.Normalize();
  Vector oParallel_1 = oV ^ oDirection1;
  Vector oParallel_2 = oV ^ oDirection2;

  if (oParallel_1.Length() < dTolerance && oParallel_2.Length() < dTolerance) {
    return true;
  } else if (oParallel_1.Length() < dTolerance) {
    oCollisionPoint = *poNode2;
    return true;
  } else if (oParallel_2.Length() < dTolerance) {
    oCollisionPoint = *poNode1;
    return true;
  } else
    return false;
}
bool ParadisCollisionServer::DoesExactCollisionPointExist(
    Point *poNode1, const unsigned int &iConstraintType1,
    const Vector &oConstraintVector1, Point *poNode2,
    const unsigned int &iConstraintType2, const Vector &oConstraintVector2,
    Point &oExactCollisionPoint) const {
  bool bDoesExist = false;
  if ((iConstraintType1 == 1) && (iConstraintType2 == 1)) {
    bDoesExist =
        GetPlanePlaneCollisionPoint(poNode1, oConstraintVector1, poNode2,
                                    oConstraintVector2, oExactCollisionPoint);
  } else if ((iConstraintType1 == 1) && (iConstraintType2 == 2)) {
    bDoesExist =
        GetLinePlaneCollisionPoint(poNode2, oConstraintVector2, poNode1,
                                   oConstraintVector1, oExactCollisionPoint);
  } else if ((iConstraintType1 == 2) && (iConstraintType2 == 1)) {
    bDoesExist =
        GetLinePlaneCollisionPoint(poNode1, oConstraintVector1, poNode2,
                                   oConstraintVector2, oExactCollisionPoint);
  }
  return bDoesExist;
}
bool ParadisCollisionServer::GetExactCollisionPoint(
    Home_t *poHome, Node_t *poNode1, Node_t *poNode2,
    Point &oCollisionPoint) const {
  oCollisionPoint.Set(0.0, 0.0, 0.0);
  if (poNode1 == NULL) {
    return false;
  }
  if (poNode2 == NULL) {
    return false;
  }
  unsigned int iConstraintType1 = 0;
  unsigned int iConstraintType2 = 0;
  Vector oConstraintVector1;
  Vector oConstraintVector2;
  real8 dX, dY, dZ;
  real8 dTol = 1E-6;
  dX = poNode1->x - poNode2->x;
  dY = poNode1->y - poNode2->y;
  dZ = poNode1->z - poNode2->z;
  Vector lineDirection, PlaneNormal;
  lineDirection.Set(dX, dY, dZ);
  ParadisSurface::GetNodeDynamicConstraint(poNode1, iConstraintType1,
                                           oConstraintVector1);
  ParadisSurface::GetNodeDynamicConstraint(poNode2, iConstraintType2,
                                           oConstraintVector2);
  Point oPoint1(poNode1->x, poNode1->y, poNode1->z);
  // get point2 relative to point1
  Point oPoint2 = GetNearestImage(poHome->param, poNode1, poNode2);

  if ((iConstraintType1 == 1) && (iConstraintType2 == 1)) {
    return GetPlanePlaneCollisionPoint(&oPoint1, oConstraintVector1, &oPoint2,
                                       oConstraintVector2, oCollisionPoint);
  } else if ((iConstraintType1 == 1) && (iConstraintType2 == 2)) {
    return GetLinePlaneCollisionPoint(&oPoint2, oConstraintVector2, &oPoint1,
                                      oConstraintVector1, oCollisionPoint);
  } else if ((iConstraintType1 == 2) && (iConstraintType2 == 1)) {
    return GetLinePlaneCollisionPoint(&oPoint1, oConstraintVector1, &oPoint2,
                                      oConstraintVector2, oCollisionPoint);
  }
  // junjie: two points on the same line
  else if ((iConstraintType1 == 2) && (iConstraintType2 == 2)) {
    return GetLineLineCollisionPoint(&oPoint1, oConstraintVector1, &oPoint2,
                                     oConstraintVector2, oCollisionPoint);
  }
  // junjie
  // junjie: enable the junctions to fully unzip
  // differentiate the unzipping from improper out of plane collisions
  else if ((iConstraintType1 == 3) && (iConstraintType2 == 1)) {
    PlaneNormal.Set(poNode2->nx[0], poNode2->ny[0], poNode2->nz[0]);
    if (fabs(lineDirection * PlaneNormal) < dTol) {
      oCollisionPoint = oPoint1;
      return true;
    }
    else return false;
  } else if ((iConstraintType2 == 3) && (iConstraintType1 == 1)) {
    PlaneNormal.Set(poNode1->nx[0], poNode1->ny[0], poNode1->nz[0]);
    if (fabs(lineDirection * PlaneNormal) < dTol) {
      oCollisionPoint = oPoint2;
      return true;
    }
    else return false;
  }
  // junjie
  return false;
}
void ParadisCollisionServer::HandleAnnihilations(Home_t *poHome) {
  unsigned int i = 0;
  Node_t *poNode1 = NULL;
  Node_t *poNode2 = NULL;
  Node_t *poNode3 = NULL;
  Node_t *poNode4 = NULL;
  bool bHasCollided =
      false; // because every node is allowed one collision a time
  int iSourceCellIndex = 0;
  int iTargetCellIndex = 0;
  int iXCell = 0;
  int iYCell = 0;
  int iZCell = 0;
  int j = 0;
  int k = 0;
  int l = 0;
  int m = 0;
  int n = 0;
  int iTargetNodeIndex = 0;
  int iDomainID = poHome->myDomain;
  double dTolerance = 1.0E-6;
  double dTemp = 0.0;
  Vector oNormal1;
  Vector oNormal2;
  Vector oBurgers1;
  Vector oBurgers2;
  Vector oSegment1;
  Vector oSegment2;

  double dMinimumLengthRatio = 0.7;
  double dCenterSpacingFactor = 0.3;
  double dLength1 = 0.0;
  double dLength2 = 0.0;
  double dMaximumCenterSpacing = 0.0;
  double dCenterSpacing = 0.0;
  double dAnnihilationDistance = poHome->param->AnnihilationDistance;
  double dMinimumAlignment = 0.9;
  double dAlignment = 0.0;
  double dNorm = 0.0;
  Point oNode1;
  Point oNode2;
  Point oNode3;
  Point oNode4;
  // open the annihlated segments tracker
  AnnihilationSegment *poSegment = NULL;
  for (i = 0; i < poHome->newNodeKeyPtr; i++) {
    poNode1 = poHome->nodeKeys[i];
    if (!IsNodeCollisionValid(poNode1)) {
      continue;
    }
    if (poNode1->numNbrs > 2) {
      continue;
    }
    bHasCollided = false;
    iSourceCellIndex = poNode1->cell2Idx;
    if (iSourceCellIndex < 0) {
      continue;
    }
    // loop over the nodes in the neighboring cells
    DecodeCell2Idx(poHome, iSourceCellIndex, &iXCell, &iYCell, &iZCell);
    for (j = iXCell - 1; j <= iXCell + 1; j++) {
      if (bHasCollided) {
        break;
      }
      for (k = iYCell - 1; k <= iYCell + 1; k++) {
        if (bHasCollided) {
          break;
        }
        for (l = iZCell - 1; l <= iZCell + 1; l++) {
          if (bHasCollided) {
            break;
          }
          iTargetCellIndex = EncodeCell2Idx(poHome, j, k, l);
          // loop over the nodes of the target cell
          iTargetNodeIndex = poHome->cell2[iTargetCellIndex];
          while (iTargetNodeIndex >= 0) {
            if (bHasCollided) {
              break;
            }
            poNode3 = poHome->cell2QentArray[iTargetNodeIndex].node;
            iTargetNodeIndex = poHome->cell2QentArray[iTargetNodeIndex].next;
            if (!IsNodeCollisionValid(poNode3)) {
              continue;
            }
            if (poNode3->numNbrs > 2) {
              continue;
            }
            if (poNode3->myTag.domainID != iDomainID) {
              continue;
            }
            if (CollisionNodeOrder(poHome, &poNode1->myTag, &poNode3->myTag) >=
                0) {
              continue;
            }
            // now the distinct arms of nodes 1 and 3 may annihilate, check for
            // that
            for (m = 0; m < poNode1->numNbrs; m++) {
              if (bHasCollided) {
                break;
              }
              poNode2 = GetNodeFromTag(poHome, poNode1->nbrTag[m]);
              if (!IsNodeCollisionValid(poNode2)) {
                continue;
              }
              if (poNode2->numNbrs > 2) {
                continue;
              }
              if (CollisionNodeOrder(poHome, &poNode1->myTag, &poNode2->myTag) >
                  0) {
                continue;
              }
              if ((poNode2->myTag.domainID == poNode3->myTag.domainID) &&
                  (poNode2->myTag.index == poNode3->myTag.index)) {
                continue;
              }
              if (!DomainOwnsSeg(poHome, OPCLASS_COLLISION, iDomainID,
                                 &poNode2->myTag)) {
                continue;
              }
              for (n = 0; n < poNode3->numNbrs; n++) {
                if (bHasCollided) {
                  break;
                }
                poNode4 = GetNodeFromTag(poHome, poNode3->nbrTag[n]);
                if (!IsNodeCollisionValid(poNode4)) {
                  continue;
                }
                if (poNode4->numNbrs > 2) {
                  continue;
                }
                if ((poNode4->myTag.domainID == poNode1->myTag.domainID) &&
                    (poNode4->myTag.index == poNode1->myTag.index)) {
                  continue;
                }
                if ((poNode4->myTag.domainID == poNode2->myTag.domainID) &&
                    (poNode4->myTag.index == poNode2->myTag.index)) {
                  continue;
                }
                if (CollisionNodeOrder(poHome, &poNode3->myTag,
                                       &poNode4->myTag) > 0) {
                  continue;
                }
                if (!DomainOwnsSeg(poHome, OPCLASS_COLLISION, iDomainID,
                                   &poNode4->myTag)) {
                  continue;
                }
                // get all the nodes in the node1 space for annihilation checks
                Point oNode1(poNode1->x, poNode1->y, poNode1->z);
                Point oNode2 = GetNearestImage(poHome->param, poNode1, poNode2);
                Point oNode3 = GetNearestImage(poHome->param, poNode1, poNode3);
                Point oNode4 = GetNearestImage(poHome->param, poNode1, poNode4);
                // now these 2 segments will annihilate if
                // 1. their lengths are roughly equal
                dLength1 = oNode1.Distance(oNode2);
                dLength2 = oNode3.Distance(oNode4);
                if ((dLength1 / dLength2 < dMinimumLengthRatio) &&
                    (dLength2 / dLength1 <
                     dMinimumLengthRatio)) // junjie: should be ||
                {
                  continue;
                }
                // 2. their centers are close enough to each other
                dMaximumCenterSpacing =
                    dCenterSpacingFactor * min(dLength1, dLength2);
                dCenterSpacing =
                    ((oNode1 + oNode2) * 0.5).Distance((oNode3 + oNode4) * 0.5);
                if (dCenterSpacing > dMaximumCenterSpacing) {
                  continue;
                }
                // 3. they on parallel planes
                oNormal1.Set(poNode1->nx[m], poNode1->ny[m], poNode1->nz[m]);
                oNormal2.Set(poNode3->nx[n], poNode3->ny[n], poNode3->nz[n]);
                if (!(oNormal1 ^ oNormal2).Length() < dTolerance) {
                  continue;
                }
                // 4. but the planes are different, and the distance between
                // them is less than the annihilation distance
                dTemp = (poNode3->x - poNode1->x) * poNode1->nx[m];
                dTemp = dTemp + (poNode3->y - poNode1->y) * poNode1->ny[m];
                dTemp = dTemp + (poNode3->z - poNode1->z) * poNode1->nz[m];
                dNorm = poNode1->nx[m] * poNode1->nx[m] +
                        poNode1->ny[m] * poNode1->ny[m] +
                        poNode1->nz[m] * poNode1->nz[m];
                dNorm = sqrt(dNorm);
                dTemp = fabs(dTemp) / dNorm;
                if (dTemp > dAnnihilationDistance) {
                  continue;
                }
                // 5. they are in the same direction up to a certain tolerance
                oSegment1.SetByPoints(oNode1, oNode2);
                oSegment2.SetByPoints(oNode3, oNode4);
                oSegment1.Normalize();
                oSegment2.Normalize();
                dAlignment = oSegment1 * oSegment2;
                if (fabs(dAlignment) < dMinimumAlignment) {
                  continue;
                }
                // 6. they have opposite Burgers vectors
                oBurgers1.Set(poNode1->burgX[m], poNode1->burgY[m],
                              poNode1->burgZ[m]);
                oBurgers2.Set(poNode3->burgX[n], poNode3->burgY[n],
                              poNode3->burgZ[n]);
                if (dAlignment > 0.0) {
                  if (!(oBurgers1 + oBurgers2).Length() < dTolerance) {
                    continue;
                  }
                } else {
                  if (!(oBurgers1 - oBurgers2).Length() < dTolerance) {
                    continue;
                  }
                }
                // if we reach this point, then the two segments can annihilate
                // simply, remove the segments connecting nodes 1-2 and 3-4
                // and pin the 4 nodes, but first, create the appropriate
                // annihilation segments and add them to the list
                poSegment = new AnnihilationSegment;
                poSegment->SetPoints(Point(poNode1->x, poNode1->y, poNode1->z),
                                     Point(poNode2->x, poNode2->y, poNode2->z));
                poSegment->SetBurgersVector(Vector(
                    poNode1->burgX[m], poNode1->burgY[m], poNode1->burgZ[m]));
                poSegment->SetTimeStep(poHome->cycle);
                m_lpoAnnihilatedSegments.push_back(poSegment);
                poSegment = new AnnihilationSegment;
                poSegment->SetPoints(Point(poNode3->x, poNode3->y, poNode3->z),
                                     Point(poNode4->x, poNode4->y, poNode4->z));
                poSegment->SetBurgersVector(Vector(
                    poNode3->burgX[n], poNode1->burgY[n], poNode1->burgZ[n]));
                poSegment->SetTimeStep(poHome->cycle);
                m_lpoAnnihilatedSegments.push_back(poSegment);

                ChangeArmBurg(poHome, poNode1, &poNode2->myTag, 0, 0, 0, 0, 0,
                              0, 1, DEL_SEG_HALF, true);
                ChangeArmBurg(poHome, poNode2, &poNode1->myTag, 0, 0, 0, 0, 0,
                              0, 1, DEL_SEG_HALF, true);
                ChangeArmBurg(poHome, poNode3, &poNode4->myTag, 0, 0, 0, 0, 0,
                              0, 1, DEL_SEG_HALF, true);
                ChangeArmBurg(poHome, poNode4, &poNode3->myTag, 0, 0, 0, 0, 0,
                              0, 1, DEL_SEG_HALF, true);
                ChangeNodeType(poHome, poNode1, PINNED_NODE, 1);
                ChangeNodeType(poHome, poNode2, PINNED_NODE, 1);
                ChangeNodeType(poHome, poNode3, PINNED_NODE, 1);
                ChangeNodeType(poHome, poNode4, PINNED_NODE, 1);
              }
            }
          }
        }
      }
    }
  }
}
void ParadisCollisionServer::HandleSegmentSegmentCollisions(
    Home_t *poHome) const {
  unsigned int i = 0;
  Node_t *poNode1 = NULL;
  Node_t *poNode2 = NULL;
  Node_t *poNode3 = NULL;
  Node_t *poNode4 = NULL;
  bool bHasCollided =
      false; // because every node is allowed one collision a time
  int iSourceCellIndex = 0;
  int iTargetCellIndex = 0;
  int iXCell = 0;
  int iYCell = 0;
  int iZCell = 0;
  int j = 0;
  int k = 0;
  int l = 0;
  int m = 0;
  int n = 0;
  int p = 0;
  int iTargetNodeIndex = 0;
  int iDomainID = poHome->myDomain;
  double dTolerance = 1.0E-6;
  double dTemp = 0.0;
  double dLocalCoordainte1 = 0.0;
  double dLocalCoordainte2 = 0.0;
  Vector oNormal1;
  Vector oNormal2;
  Point oCollisionPoint;
  Point oNearPoint1;
  Point oNearPoint2;
  double dNewPosition[3];
  Node_t *poMergeNode1 = NULL;
  Node_t *poMergeNode2 = NULL;
  Node_t *poMergeNeighbour = NULL;
  Node_t *poTargetNode = NULL;
  int iMergeStatus = 0;
  Vector oBurgers1;
  Vector oBurgers2;
  Vector oSegment1;
  Vector oSegment2;
  double dCloseNodeThreshold = 1.0;
  bool bExactCollisionPointExists = false;
  bool bNo1Split = false;
  bool bNo2Split = false;

  for (i = 0; i < poHome->newNodeKeyPtr; i++) {
    poNode1 = poHome->nodeKeys[i];
    if (!IsNodeCollisionValid(poNode1)) {
      continue;
    }
    if (poNode1->numNbrs > 2) {
      if (poHome->param->EnableTwinPlaneCrossSlip == 1 &&
          poNode1->csState == AT_INTERSECTION)
        ;
      else
        continue;
    }
    bHasCollided = false;
    iSourceCellIndex = poNode1->cell2Idx;
    if (iSourceCellIndex < 0) {
      continue;
    }
    // loop over the nodes in the neighboring cells
    DecodeCell2Idx(poHome, iSourceCellIndex, &iXCell, &iYCell, &iZCell);
    for (j = iXCell - 1; j <= iXCell + 1; j++) {
      if (bHasCollided) {
        break;
      }
      for (k = iYCell - 1; k <= iYCell + 1; k++) {
        if (bHasCollided) {
          break;
        }
        for (l = iZCell - 1; l <= iZCell + 1; l++) {
          if (bHasCollided) {
            break;
          }
          iTargetCellIndex = EncodeCell2Idx(poHome, j, k, l);
          // loop over the nodes of the target cell
          iTargetNodeIndex = poHome->cell2[iTargetCellIndex];
          while (iTargetNodeIndex >= 0) {
            if (bHasCollided) {
              break;
            }
            poNode3 = poHome->cell2QentArray[iTargetNodeIndex].node;
            iTargetNodeIndex = poHome->cell2QentArray[iTargetNodeIndex].next;
            if (!IsNodeCollisionValid(poNode3)) {
              continue;
            }
            if (poNode3->numNbrs > 2) {
              if (poHome->param->EnableTwinPlaneCrossSlip == 1 &&
                  poNode3->csState == AT_INTERSECTION)
                ;
              else
                continue;
            }
            if (poNode3->myTag.domainID != iDomainID) {
              continue;
            }
            if (CollisionNodeOrder(poHome, &poNode1->myTag, &poNode3->myTag) >=
                0) {
              continue;
            }
            // now the distinct arms of nodes 1 and 3 may collide, check for
            // that
            for (m = 0; m < poNode1->numNbrs; m++) {
              if (bHasCollided) {
                break;
              }
              poNode2 = GetNodeFromTag(poHome, poNode1->nbrTag[m]);
              if (!IsNodeCollisionValid(poNode2)) {
                continue;
              }
              if (poNode2->numNbrs > 2) {
                if (poHome->param->EnableTwinPlaneCrossSlip == 1 &&
                    poNode2->csState == AT_INTERSECTION)
                  ;
                else
                  continue;
              }
              if (CollisionNodeOrder(poHome, &poNode1->myTag, &poNode2->myTag) >
                  0) {
                continue;
              }
              if ((poNode2->myTag.domainID == poNode3->myTag.domainID) &&
                  (poNode2->myTag.index == poNode3->myTag.index)) {
                continue;
              }
              if (!DomainOwnsSeg(poHome, OPCLASS_COLLISION, iDomainID,
                                 &poNode2->myTag)) {
                continue;
              }
              for (n = 0; n < poNode3->numNbrs; n++) {
                if (bHasCollided) {
                  break;
                }
                poNode4 = GetNodeFromTag(poHome, poNode3->nbrTag[n]);
                if (!IsNodeCollisionValid(poNode4)) {
                  continue;
                }
                if (poNode4->numNbrs > 2) {
                  if (poHome->param->EnableTwinPlaneCrossSlip == 1 &&
                      poNode4->csState == AT_INTERSECTION)
                    ;
                  else
                    continue;
                }
                if ((poNode4->myTag.domainID == poNode1->myTag.domainID) &&
                    (poNode4->myTag.index == poNode1->myTag.index)) {
                  continue;
                }
                if ((poNode4->myTag.domainID == poNode2->myTag.domainID) &&
                    (poNode4->myTag.index == poNode2->myTag.index)) {
                  continue;
                }
                if (CollisionNodeOrder(poHome, &poNode3->myTag,
                                       &poNode4->myTag) > 0) {
                  continue;
                }
                if (!DomainOwnsSeg(poHome, OPCLASS_COLLISION, iDomainID,
                                   &poNode4->myTag)) {
                  continue;
                }
                // are segments close enough for collision ?? are they
                // approaching each other or moving away from each other ??
                if (!AreSegmentsCloseAndApproaching(poHome, poNode1, poNode2,
                                                    poNode3, poNode4,
                                                    oNearPoint1, oNearPoint2)) {
                  continue;
                }
                if (poHome->param->EnableTwinPlaneCrossSlip == 1) {
                  // junjie:only allow collision when:
                  // 1. neither segments are on the twin plane (all the 4 nodes
                  // involved are off the twin plane)
                  // 2. both segments are on the twin plane (all the 4 nodes
                  // involved are on the twin plane)
                  double A, B, C, D;
                  double X_1, Y_1, Z_1;
                  double dTol =
                      1E-8; // set a tolerance for geometrical detection
                  Point pNode1, pNode2, pNode3, pNode4;
                  pNode1.Set(poNode1->x, poNode1->y, poNode1->z);
                  pNode2.Set(poNode2->x, poNode2->y, poNode2->z);
                  pNode3.Set(poNode3->x, poNode3->y, poNode3->z);
                  pNode4.Set(poNode4->x, poNode4->y, poNode4->z);
                  A = poHome->param->A;
                  B = poHome->param->B;
                  C = poHome->param->C;
                  X_1 = poHome->param->X_1;
                  Y_1 = poHome->param->Y_1;
                  Z_1 = poHome->param->Z_1;
                  ;
                  D = -(A * X_1 + B * Y_1 + C * Z_1);
                  Plane crossSlipPlane;
                  crossSlipPlane.Set(Vector(A, B, C), Point(X_1, Y_1, Z_1));
                  if ((crossSlipPlane.GetPointDistance(pNode1) < dTol) ||
                      (crossSlipPlane.GetPointDistance(pNode2) < dTol) ||
                      (crossSlipPlane.GetPointDistance(pNode3) < dTol) ||
                      (crossSlipPlane.GetPointDistance(pNode4) < dTol)) {
                    if ((crossSlipPlane.GetPointDistance(pNode1) >= dTol) ||
                        (crossSlipPlane.GetPointDistance(pNode2) >= dTol) ||
                        (crossSlipPlane.GetPointDistance(pNode3) >= dTol) ||
                        (crossSlipPlane.GetPointDistance(pNode4) >= dTol)) {
                      continue;
                    }
                  }
                  // if
                  // (!SegSegCanHappen(poHome,poNode1,poNode2,poNode3,poNode4))
                  // {
                  //	continue;
                  //}
                  // junjie
                }

                // now these 2 segments are supposed to collide, make sure that
                // the extra kinematic constraints will not cause problems
                oNormal1.Set(poNode1->nx[m], poNode1->ny[m], poNode1->nz[m]);
                oNormal2.Set(poNode3->nx[n], poNode3->ny[n], poNode3->nz[n]);
                oNormal1.Normalize();
                oNormal2.Normalize();

                // now the segments can definitely collide, but if the nearest
                // point on either of the two segments is too close to one of
                // the segment nodes, skip the collision this time, because we
                // don't want to split the nodes and place the new nodes at the
                // same position. before trying anything, just make sure that an
                // exact collision point exists for the nodes to be created the
                // constraint type will always be 1 in this case because this
                // check will be of use only if the nodes are the newly split
                // nodes
                bExactCollisionPointExists = DoesExactCollisionPointExist(
                    &oNearPoint1, 1, oNormal1, &oNearPoint2, 1, oNormal2,
                    oCollisionPoint);
                // if the collision node is too far from either nodes, skip the
                // collision, use rmin to be the allowed collision distance
                if (oNearPoint1.Distance(oCollisionPoint) >
                    poHome->param->rmax) {
                  continue;
                }
                if (oNearPoint2.Distance(oCollisionPoint) >
                    poHome->param->rmax) {
                  continue;
                }

                bNo1Split = false;
                bNo2Split = false;
                if (GetNearestImage(poHome->param, poNode1, oNearPoint1)
                        .Distance(Point(poNode1->x, poNode1->y, poNode1->z)) <
                    dCloseNodeThreshold) {
                  poMergeNode1 = poNode1;
                  bNo1Split = true;
                } else if (GetNearestImage(poHome->param, poNode2, oNearPoint1)
                               .Distance(
                                   Point(poNode2->x, poNode2->y, poNode2->z)) <
                           dCloseNodeThreshold) {
                  poMergeNode1 = poNode2;
                  bNo1Split = true;
                } else if (bExactCollisionPointExists) {
                  // split the first arm
                  dNewPosition[0] = oNearPoint1.GetX();
                  dNewPosition[1] = oNearPoint1.GetY();
                  dNewPosition[2] = oNearPoint1.GetZ();
                  SplitNode(poHome, poNode1, dNewPosition, m, 1, poMergeNode1);
                } else {
                  continue;
                }

                if (GetNearestImage(poHome->param, poNode3, oNearPoint2)
                        .Distance(Point(poNode3->x, poNode3->y, poNode3->z)) <
                    dCloseNodeThreshold) {
                  poMergeNode2 = poNode3;
                  bNo2Split = true;
                } else if (GetNearestImage(poHome->param, poNode4, oNearPoint2)
                               .Distance(
                                   Point(poNode4->x, poNode4->y, poNode4->z)) <
                           dCloseNodeThreshold) {
                  poMergeNode2 = poNode4;
                  bNo2Split = true;
                } else if (bExactCollisionPointExists) {
                  // split the second arm
                  dNewPosition[0] = oNearPoint2.GetX();
                  dNewPosition[1] = oNearPoint2.GetY();
                  dNewPosition[2] = oNearPoint2.GetZ();
                  SplitNode(poHome, poNode3, dNewPosition, n, 1, poMergeNode2);
                } else {
                  // this condition should be unreachable
                  continue;
                }

                if (bNo1Split || bNo2Split) {
                  if (!GetExactCollisionPoint(poHome, poMergeNode1,
                                              poMergeNode2, oCollisionPoint)) {
                    continue;
                  }
                }
                dNewPosition[0] = oCollisionPoint.GetX();
                dNewPosition[1] = oCollisionPoint.GetY();
                dNewPosition[2] = oCollisionPoint.GetZ();
                iMergeStatus =
                    MergeNode(poHome, OPCLASS_COLLISION, poMergeNode1,
                              poMergeNode2, dNewPosition, &poTargetNode, 1);
                if ((iMergeStatus & MERGE_SUCCESS) == 0) {
                  continue;
                }
                if (poTargetNode == NULL) {
                  continue;
                }
                poTargetNode->flags |= NO_COLLISIONS;
                bHasCollided = true;
              }
            }
          }
        }
      }
    }
  }
  fflush(NULL);
}
void ParadisCollisionServer::HandleHingeJointCollisions(Home_t *poHome) const {
  unsigned int i = 0;
  unsigned int j = 0;
  unsigned int k = 0;
  unsigned int l = 0;
  Node_t *poNode1 = NULL;
  Node_t *poNode2 = NULL;
  Node_t *poNode3 = NULL;
  Node_t *poCheckNode = NULL;
  Node_t *poTargetNode = NULL;
  double dNewPosition[3];
  Point oCollisionPoint;
  Node_t *poMergeNode1 = NULL;
  Node_t *poMergeNode2 = NULL;
  int iSplitArmNode = 0;
  int iMergeStatus = 0;
  bool bTerminalsConnected = false;

  for (i = 0; i < poHome->newNodeKeyPtr; i++) {
    poNode1 = poHome->nodeKeys[i];
    if (!IsNodeCollisionValid(poNode1)) {
      continue;
    }
    for (j = 0; j < poNode1->numNbrs; j++) {
      if (poNode1->myTag.domainID != poNode1->nbrTag[j].domainID) {
        continue;
      }
      for (k = j + 1; k < poNode1->numNbrs; k++) {
        if (poNode1->myTag.domainID != poNode1->nbrTag[k].domainID) {
          continue;
        }
        poNode2 = GetNodeFromTag(poHome, poNode1->nbrTag[j]);
        if (!IsNodeCollisionValid(poNode2)) {
          continue;
        }
        poNode3 = GetNodeFromTag(poHome, poNode1->nbrTag[k]);
        if (!IsNodeCollisionValid(poNode3)) {
          continue;
        }
        // make sure of 2 things
        // 1. the nodes 2 and 3 are not the same (they cannot be anyway, but
        // just to be on the safe side
        if (poNode2->myTag.index == poNode3->myTag.index) {
          continue;
        }
        // 2. the nodes 2 and 3 are not connected
        bTerminalsConnected = false;
        for (l = 0; l < poNode2->numNbrs; l++) {
          poCheckNode = GetNodeFromTag(poHome, poNode2->nbrTag[l]);
          if (poCheckNode == NULL) {
            continue;
          }
          if (poCheckNode->myTag.index == poNode3->myTag.index) {
            if (poCheckNode->myTag.domainID == poNode3->myTag.domainID) {
              bTerminalsConnected = true;
              break;
            }
          }
        }
        if (bTerminalsConnected) {
          continue;
        }
        // now we can make the further kinematic checks
        // junjie
        // if (!HingeJointCanHappen(poHome,poNode1,poNode2,poNode3)) {
        //	continue;
        //}
        // junjie
        if (!IsTriangleCollapsing(poHome, poNode1, poNode2, poNode3)) {
          continue;
        }
        if (poHome->param->EnableTwinPlaneCrossSlip == 1) {
          // junjie:only allow collision when:
          // poNode2 & poNode3 are both on or off the csplane
          double A, B, C, D;
          double X_1, Y_1, Z_1;
          double dTol = 1E-8; // set a tolerance for geometrical detection
          Point pNode1, pNode2, pNode3;
          pNode1.Set(poNode1->x, poNode1->y, poNode1->z);
          pNode2.Set(poNode2->x, poNode2->y, poNode2->z);
          pNode3.Set(poNode3->x, poNode3->y, poNode3->z);
          A = poHome->param->A;
          B = poHome->param->B;
          C = poHome->param->C;
          X_1 = poHome->param->X_1;
          Y_1 = poHome->param->Y_1;
          Z_1 = poHome->param->Z_1;
          ;
          D = -(A * X_1 + B * Y_1 + C * Z_1);
          Plane crossSlipPlane;
          crossSlipPlane.Set(Vector(A, B, C), Point(X_1, Y_1, Z_1));
          if ((crossSlipPlane.GetPointDistance(pNode2) >= dTol &&
               crossSlipPlane.GetPointDistance(pNode3) >= dTol) ||
              (crossSlipPlane.GetPointDistance(pNode2) < dTol &&
               crossSlipPlane.GetPointDistance(pNode3) < dTol))
            ;
          else
            continue;
          // junjie
        }

        // merge the two end nodes if possible, otherwise, skip the collision
        poMergeNode1 = poNode2;
        poMergeNode2 = poNode3;
        if (!GetExactCollisionPoint(poHome, poMergeNode1, poMergeNode2,
                                    oCollisionPoint)) {
          continue;
        }
        dNewPosition[0] = oCollisionPoint.GetX();
        dNewPosition[1] = oCollisionPoint.GetY();
        dNewPosition[2] = oCollisionPoint.GetZ();

        iMergeStatus = MergeNode(poHome, OPCLASS_COLLISION, poMergeNode1,
                                 poMergeNode2, dNewPosition, &poTargetNode, 1);

        if ((iMergeStatus & MERGE_SUCCESS) == 0) {
          continue;
        }
        // now since the merge worked, make sure that any apb points were
        // properly updated
        if (poTargetNode != NULL) {
          poTargetNode->flags |= NO_COLLISIONS;
        }
      }
    }
  }
}
void ParadisCollisionServer::HandleTriangularLoops(Home_t *poHome) const {
  unsigned int i = 0;
  unsigned int j = 0;
  unsigned int k = 0;
  unsigned int l = 0;
  Node_t *poNode1 = NULL;
  Node_t *poNode2 = NULL;
  Node_t *poNode3 = NULL;
  Node_t *poMergeNode1 = NULL;
  Node_t *poMergeNode2 = NULL;
  Node_t *poTargetNode = NULL;
  bool bAreNeighbours = false;
  bool bFoundMergeNodes = false;
  double dNewPosition[3];
  int iMergeStatus = 0;
  for (i = 0; i < poHome->newNodeKeyPtr; i++) {
    poNode1 = poHome->nodeKeys[i];
    if (!IsNodeCollisionValid(poNode1)) {
      continue;
    }
    for (j = 0; j < poNode1->numNbrs; j++) {
      if (poNode1->myTag.domainID != poNode1->nbrTag[j].domainID) {
        continue;
      }
      for (k = j + 1; k < poNode1->numNbrs; k++) {
        if (poNode1->myTag.domainID != poNode1->nbrTag[k].domainID) {
          continue;
        }
        poNode2 = GetNodeFromTag(poHome, poNode1->nbrTag[j]);
        if (!IsNodeCollisionValid(poNode2)) {
          continue;
        }
        poNode3 = GetNodeFromTag(poHome, poNode1->nbrTag[k]);
        if (!IsNodeCollisionValid(poNode3)) {
          continue;
        }
        // check to see if the nodes 2 and 3 are neighbours
        bAreNeighbours = false;
        for (l = 0; l < poNode2->numNbrs; l++) {
          if (poNode2->nbrTag[l].index == poNode3->myTag.index) {
            if (poNode2->nbrTag[l].domainID == poNode3->myTag.domainID) {
              bAreNeighbours = true;
              break;
            }
          }
        }
        if (!bAreNeighbours) {
          continue;
        }
        // now we have a triangle, make sure that one of the 3 nodes has exactly
        // 2 arms this node is to be removed
        poMergeNode1 = NULL;
        poMergeNode2 = NULL;
        bFoundMergeNodes = false;
        if (poNode1->numNbrs == 2) {
          poMergeNode1 = poNode1;
          if ((poNode2->flags & NO_MESH_COARSEN) == 0) {
            poMergeNode2 = poNode2;
            bFoundMergeNodes = true;
          } else if ((poNode3->flags & NO_MESH_COARSEN) == 0) {
            poMergeNode2 = poNode3;
            bFoundMergeNodes = true;
          }
        }
        if (!bFoundMergeNodes) {
          if (poNode2->numNbrs == 2) {
            poMergeNode1 = poNode2;
            if ((poNode1->flags & NO_MESH_COARSEN) == 0) {
              poMergeNode2 = poNode1;
              bFoundMergeNodes = true;
            } else if ((poNode3->flags & NO_MESH_COARSEN) == 0) {
              poMergeNode2 = poNode3;
              bFoundMergeNodes = true;
            }
          }
        }
        if (!bFoundMergeNodes) {
          if (poNode3->numNbrs == 2) {
            poMergeNode1 = poNode3;
            if ((poNode1->flags & NO_MESH_COARSEN) == 0) {
              poMergeNode2 = poNode1;
              bFoundMergeNodes = true;
            } else if ((poNode2->flags & NO_MESH_COARSEN) == 0) {
              poMergeNode2 = poNode2;
              bFoundMergeNodes = true;
            }
          }
        }
        if (!bFoundMergeNodes) {
          continue;
        }
        // now everything is set, we have the two nodes to be merged, merge and
        // place in the position of the second merge node
        dNewPosition[0] = poMergeNode2->x;
        dNewPosition[1] = poMergeNode2->y;
        dNewPosition[2] = poMergeNode2->z;

        iMergeStatus = MergeNode(poHome, OPCLASS_COLLISION, poMergeNode1,
                                 poMergeNode2, dNewPosition, &poTargetNode, 1);
        if ((iMergeStatus & MERGE_SUCCESS) == 0) {
          continue;
        }
        if (poTargetNode != NULL) {
          poTargetNode->flags |= NO_COLLISIONS;
          poTargetNode->flags |= NO_MESH_COARSEN;
        }
      }
    }
  }
}
void ParadisCollisionServer::WriteAnnihilatedSegments() {
  // write any segments in the list and delete them as we write
  list<AnnihilationSegment *>::iterator liSegments;
  for (liSegments = m_lpoAnnihilatedSegments.begin();
       liSegments != m_lpoAnnihilatedSegments.end(); liSegments++) {
    if ((*liSegments) != NULL) {
      (*liSegments)->Write(m_fpAnnihilatedSegments);
      delete (*liSegments);
    }
  }
  m_lpoAnnihilatedSegments.clear();
}
